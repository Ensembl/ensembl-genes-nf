#!/usr/bin/env python3
"""
Merge study matrices into a global matrix with Zarr output using streaming merge.

Takes all study matrix outputs (optionally for a single partition) and produces:
1. Global Zarr matrix (chunked by reads for efficient locus queries)
2. Global unique reads FASTA (for single-pass alignment)
3. Metadata parquet dataset (read_id, length, etc.)

Uses k-way sequence streaming for O(n log k) complexity.
When used with partitioned inputs, each invocation handles one dinucleotide partition.
Designed for 100-250 studies, ~100-200M total unique reads.
"""

import argparse
import contextlib
import gzip
import heapq
import json
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

import numpy as np
import scipy.sparse as sp

try:
    import zarr
    import numcodecs
    HAVE_ZARR = True
except ImportError:
    HAVE_ZARR = False
    print("Warning: zarr not available, will output NPZ instead", file=sys.stderr)

try:
    import polars as pl
    HAVE_POLARS = True
except ImportError:
    HAVE_POLARS = False

try:
    import xxhash
    USE_XXHASH = True
except ImportError:
    import hashlib
    USE_XXHASH = False


def hash_sequence(seq: str) -> int:
    """Hash a sequence to a 64-bit integer."""
    if USE_XXHASH:
        return xxhash.xxh64(seq.encode()).intdigest()
    else:
        return int(hashlib.md5(seq.encode()).hexdigest()[:16], 16)


VocabEntry = Union[int, List[Tuple[str, int]]]


@dataclass
class StudyInput:
    study_id: str
    sequences_path: Path
    matrix_path: Path
    metadata_path: Path
    n_reads: int
    samples: List[str]


class MetadataShardWriter:
    """Write read metadata as bounded TSV shards, then convert to parquet."""

    def __init__(self, output_path: Path, temp_dir: Path, shard_rows: int):
        self.output_path = Path(output_path)
        self.temp_dir = Path(temp_dir)
        self.shard_rows = max(1, shard_rows)
        self.shard_index = 0
        self.rows_in_shard = 0
        self.total_rows = 0
        self.current_handle = None
        self.current_path: Optional[Path] = None
        self.parquet_parts: List[Path] = []

    def __enter__(self):
        self._open_next_shard()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def _open_next_shard(self):
        self._convert_current_shard()
        shard_name = f"metadata_part_{self.shard_index:05d}"
        self.current_path = self.temp_dir / f"{shard_name}.tsv"
        self.current_handle = open(self.current_path, 'w')
        self.current_handle.write("read_id\tlength\n")
        self.rows_in_shard = 0
        self.shard_index += 1

    def _close_handle(self):
        self.close()

    def _convert_current_shard(self):
        if self.current_path is None:
            return

        self._close_handle()
        if self.rows_in_shard == 0:
            try:
                self.current_path.unlink()
            except FileNotFoundError:
                pass
            self.current_path = None
            return

        if not HAVE_POLARS:
            raise RuntimeError("polars is required to stream global metadata parquet")

        parquet_path = self.temp_dir / f"metadata_part_{len(self.parquet_parts):05d}.parquet"
        self._scan_tsv(self.current_path).sink_parquet(
            parquet_path,
            compression='zstd',
        )
        try:
            self.current_path.unlink()
        except FileNotFoundError:
            pass
        self.parquet_parts.append(parquet_path)
        self.current_path = None
        self.rows_in_shard = 0

    def write(self, read_id: int, length: int):
        if self.current_handle is None:
            self._open_next_shard()
        if self.rows_in_shard >= self.shard_rows:
            self._open_next_shard()

        self.current_handle.write(f"{read_id}\t{length}\n")
        self.rows_in_shard += 1
        self.total_rows += 1

    def close(self):
        if self.current_handle is not None:
            self.current_handle.close()
            self.current_handle = None

    def _scan_tsv(self, tsv_path: Path):
        scan_kwargs = {'separator': '\t'}
        try:
            return pl.scan_csv(
                tsv_path,
                **scan_kwargs,
                schema_overrides={'read_id': pl.UInt64, 'length': pl.UInt32},
            )
        except TypeError:
            return pl.scan_csv(
                tsv_path,
                **scan_kwargs,
                dtypes={'read_id': pl.UInt64, 'length': pl.UInt32},
            )

    def finalize(self) -> str:
        """Convert shards to parquet. Returns 'file' or 'dataset'."""
        self._convert_current_shard()
        if not HAVE_POLARS:
            raise RuntimeError("polars is required to stream global metadata parquet")

        if self.output_path.exists():
            if self.output_path.is_dir():
                shutil.rmtree(self.output_path)
            else:
                self.output_path.unlink()

        if self.total_rows == 0:
            pl.DataFrame(
                {
                    'read_id': pl.Series([], dtype=pl.UInt64),
                    'length': pl.Series([], dtype=pl.UInt32),
                }
            ).write_parquet(self.output_path)
            return 'file'

        if len(self.parquet_parts) == 1:
            shutil.move(str(self.parquet_parts[0]), self.output_path)
            metadata_format = 'file'
        else:
            self.output_path.mkdir(parents=True, exist_ok=True)
            for idx, parquet_path in enumerate(self.parquet_parts):
                shutil.move(str(parquet_path), self.output_path / f"part-{idx:05d}.parquet")
            metadata_format = 'dataset'

        return metadata_format


def resolve_local_id(seq_to_id: Dict[int, VocabEntry], seq: str) -> Optional[int]:
    """
    Resolve a study-local row ID for a sequence.

    New vocab files store collision buckets when multiple sequences share the
    same 64-bit hash. Older vocab files store hash -> int only; those still load,
    but cannot disambiguate a hash collision that was not recorded at build time.
    """
    entry = seq_to_id.get(hash_sequence(seq))
    if entry is None:
        return None
    if isinstance(entry, int):
        return entry

    for candidate_seq, local_id in entry:
        if candidate_seq == seq:
            return local_id
    return None


# =============================================================================
# Streaming merge helpers
# =============================================================================

def read_study_sequences(sequences_path: Path, study_id: str) -> Iterator[Tuple[str, int]]:
    """
    Stream sorted sequences from a study's gzipped sequences file.
    Yields: (sequence, local_id)
    """
    previous_seq: Optional[str] = None
    local_id = 0

    with gzip.open(sequences_path, 'rt') as f:
        for line_num, line in enumerate(f, start=1):
            seq = line.strip()
            if not seq:
                raise ValueError(f"Empty sequence in {sequences_path}:{line_num}")
            if previous_seq is not None and seq <= previous_seq:
                raise ValueError(
                    f"Study sequence file is not strictly sorted for {study_id} "
                    f"at {sequences_path}:{line_num}: {seq!r} follows {previous_seq!r}"
                )

            yield seq, local_id
            previous_seq = seq
            local_id += 1


def read_intermediate_file(path: Path) -> Iterator[Tuple[str, List[Tuple[str, int]]]]:
    """
    Read intermediate merge file.
    Format: sequence\tstudy1:local_id,study2:local_id,...

    Yields: (sequence, [(study_id, local_id), ...])
    """
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            seq = parts[0]
            occurrences = []
            if len(parts) > 1 and parts[1]:
                for item in parts[1].split(','):
                    study_id, local_id = item.rsplit(':', 1)
                    occurrences.append((study_id, int(local_id)))
            yield seq, occurrences


def write_intermediate_file(path: Path, data: Iterator[Tuple[str, List[Tuple[str, int]]]]):
    """
    Write intermediate merge file.
    Format: sequence\tstudy1:local_id,study2:local_id,...
    """
    with open(path, 'w') as f:
        for seq, occurrences in data:
            occurrences_str = ','.join(
                f"{study_id}:{local_id}" for study_id, local_id in occurrences
            )
            f.write(f"{seq}\t{occurrences_str}\n")


def merge_two_sources(
    source1: Iterator[Tuple[str, List[Tuple[str, int]]]],
    source2: Iterator[Tuple[str, List[Tuple[str, int]]]]
) -> Iterator[Tuple[str, List[Tuple[str, int]]]]:
    """
    Two-pointer merge of two sorted sources.
    Each source yields (sequence, [study_ids])

    This is the core merge operation - O(n) for two sorted lists.
    """
    item1: Optional[Tuple[str, List[Tuple[str, int]]]] = None
    item2: Optional[Tuple[str, List[Tuple[str, int]]]] = None

    # Get first items
    try:
        item1 = next(source1)
    except StopIteration:
        item1 = None

    try:
        item2 = next(source2)
    except StopIteration:
        item2 = None

    while item1 is not None or item2 is not None:
        if item1 is None:
            # Source 1 exhausted, drain source 2
            yield item2
            for remaining in source2:
                yield remaining
            break
        elif item2 is None:
            # Source 2 exhausted, drain source 1
            yield item1
            for remaining in source1:
                yield remaining
            break
        else:
            # Compare sequences
            seq1, occurrences1 = item1
            seq2, occurrences2 = item2

            if seq1 < seq2:
                yield item1
                try:
                    item1 = next(source1)
                except StopIteration:
                    item1 = None
            elif seq1 > seq2:
                yield item2
                try:
                    item2 = next(source2)
                except StopIteration:
                    item2 = None
            else:
                # Same sequence - merge occurrence lists
                merged_occurrences = occurrences1 + occurrences2
                yield (seq1, merged_occurrences)
                try:
                    item1 = next(source1)
                except StopIteration:
                    item1 = None
                try:
                    item2 = next(source2)
                except StopIteration:
                    item2 = None


def study_to_intermediate(sequences_path: Path, study_id: str) -> Iterator[Tuple[str, List[Tuple[str, int]]]]:
    """Convert a study's sequences to intermediate format."""
    for seq, local_id in read_study_sequences(sequences_path, study_id):
        yield (seq, [(study_id, local_id)])


class TournamentMerger:
    """
    Tournament-style pairwise merger for sorted sequence files.

    Merges N files in O(log N) rounds, each round doing pairwise merges.
    Total complexity: O(n log k) where n = total entries, k = number of files.
    """

    def __init__(self, temp_dir: Path):
        self.temp_dir = Path(temp_dir)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.temp_files: List[Path] = []
        self.round_num = 0

    def _new_temp_file(self, prefix: str = "merge") -> Path:
        """Create a new temp file path."""
        path = self.temp_dir / f"{prefix}_r{self.round_num}_{len(self.temp_files)}.txt"
        self.temp_files.append(path)
        return path

    def merge_studies(self, study_inputs: List[Tuple[str, Path]]) -> Path:
        """
        Merge multiple (study_id, sequences_path) pairs using tournament merge.

        Args:
            study_inputs: List of (study_id, sequences_gz_path) tuples

        Returns:
            Path to final merged intermediate file
        """
        if not study_inputs:
            raise ValueError("No inputs to merge")

        # Convert all inputs to intermediate format files
        current_files: List[Path] = []

        print(f"  Round 0: Converting {len(study_inputs)} study sequence files", file=sys.stderr)

        for study_id, seq_path in study_inputs:
            out_path = self._new_temp_file("init")
            write_intermediate_file(out_path, study_to_intermediate(seq_path, study_id))
            current_files.append(out_path)

        # Tournament rounds
        self.round_num = 1
        while len(current_files) > 1:
            next_files: List[Path] = []
            n_pairs = len(current_files) // 2

            print(f"  Round {self.round_num}: Merging {len(current_files)} files into {n_pairs + len(current_files) % 2}", file=sys.stderr)

            # Pairwise merges
            for i in range(0, len(current_files) - 1, 2):
                file1 = current_files[i]
                file2 = current_files[i + 1]

                out_path = self._new_temp_file("merge")

                # Stream merge
                source1 = read_intermediate_file(file1)
                source2 = read_intermediate_file(file2)
                merged = merge_two_sources(source1, source2)
                write_intermediate_file(out_path, merged)

                next_files.append(out_path)

            # Handle odd file (pass through to next round)
            if len(current_files) % 2 == 1:
                next_files.append(current_files[-1])

            current_files = next_files
            self.round_num += 1

        return current_files[0]

    def cleanup(self, keep_final: bool = True):
        """Remove temporary files."""
        for path in self.temp_files[:-1] if keep_final else self.temp_files:
            try:
                path.unlink()
            except FileNotFoundError:
                pass


# =============================================================================
# Global matrix builder
# =============================================================================

class GlobalMatrixMerger:
    """
    Merges study matrices into a global Zarr matrix using tournament merge.

    Process:
    1. Tournament merge study sequences to build global vocabulary
    2. Build global->study ID remapping
    3. Stream study matrices into global Zarr array
    """

    def __init__(
        self,
        output_dir: Path,
        chunk_size: int = 10000,
        partition: str = None,
        metadata_shard_rows: int = 100_000_000,
        write_fasta: bool = True,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.chunk_size = chunk_size
        self.partition = partition  # e.g., "AA", "AC", etc. or None
        self.metadata_shard_rows = metadata_shard_rows
        self.metadata_format = 'file'
        self.write_fasta = write_fasta

        # Output file prefix
        self.prefix = f"global.{partition}" if partition else "global"

        # Per-study info
        self.study_remaps: Dict[str, np.ndarray] = {}  # study -> disk-backed local_id -> global_id
        self.study_n_reads: Dict[str, int] = {}
        self.study_remap_counts: Dict[str, int] = {}
        self.all_samples: List[str] = []
        self.study_sample_offsets: Dict[str, int] = {}  # study -> sample column offset
        self.studies: List[str] = []
        self.n_reads = 0

        self.fasta_path = self.output_dir / f'{self.prefix}_reads.fasta'
        self.metadata_path = self.output_dir / f'{self.prefix}_metadata.parquet'

    def _prepare_study_remaps(self, study_inputs: List[StudyInput], temp_dir: Path):
        """Create disk-backed remap arrays and record sample offsets."""
        for study_idx, study in enumerate(study_inputs):
            duplicate_samples = set(study.samples).intersection(self.all_samples)
            if duplicate_samples:
                examples = ', '.join(sorted(duplicate_samples)[:5])
                raise ValueError(
                    f"Duplicate sample IDs across study matrices before {study.study_id}: {examples}"
                )

            self.study_sample_offsets[study.study_id] = len(self.all_samples)
            self.all_samples.extend(study.samples)
            self.studies.append(study.study_id)
            self.study_n_reads[study.study_id] = study.n_reads
            self.study_remap_counts[study.study_id] = 0

            if study.n_reads == 0:
                self.study_remaps[study.study_id] = np.asarray([], dtype=np.uint64)
                continue

            remap_path = temp_dir / f"study_{study_idx}.local_to_global.u64"
            self.study_remaps[study.study_id] = np.memmap(
                remap_path,
                dtype=np.uint64,
                mode='w+',
                shape=(study.n_reads,)
            )

    def _stream_global_vocabulary(self, study_inputs: List[StudyInput], temp_dir: Path):
        """
        K-way merge study sequence files and write outputs/remaps as a stream.

        This replaces the old all-vocab load. Memory now scales with number of
        studies, not number of study-local or global unique sequences.
        """
        heap: List[Tuple[str, int, int]] = []  # sequence, study_index, local_id
        iterators: List[Iterator[Tuple[str, int]]] = []

        def advance(study_index: int):
            try:
                seq, local_id = next(iterators[study_index])
            except StopIteration:
                return
            heapq.heappush(heap, (seq, study_index, local_id))

        with contextlib.ExitStack() as stack:
            if self.write_fasta:
                fasta = stack.enter_context(open(self.fasta_path, 'w'))
            else:
                self.fasta_path.write_text("")
                fasta = None
            metadata_writer = stack.enter_context(MetadataShardWriter(
                self.metadata_path,
                temp_dir,
                self.metadata_shard_rows,
            ))

            for study in study_inputs:
                iterators.append(read_study_sequences(study.sequences_path, study.study_id))
            for study_index in range(len(study_inputs)):
                advance(study_index)

            global_id = 0
            while heap:
                seq = heap[0][0]
                occurrences: List[Tuple[int, int]] = []

                while heap and heap[0][0] == seq:
                    _, study_index, local_id = heapq.heappop(heap)
                    occurrences.append((study_index, local_id))
                    advance(study_index)

                seen_studies = set()
                for study_index, local_id in occurrences:
                    study = study_inputs[study_index]
                    if study.study_id in seen_studies:
                        raise ValueError(
                            f"Duplicate sequence {seq!r} within study {study.study_id}; "
                            "study sequence files must contain one row per unique sequence"
                        )
                    seen_studies.add(study.study_id)

                    if local_id >= study.n_reads:
                        raise ValueError(
                            f"Local row {local_id} for {study.study_id} exceeds "
                            f"metadata n_reads={study.n_reads}"
                        )

                    self.study_remaps[study.study_id][local_id] = global_id
                    self.study_remap_counts[study.study_id] += 1

                if fasta is not None:
                    fasta.write(f">read_{global_id}\n{seq}\n")
                metadata_writer.write(global_id, len(seq))
                global_id += 1

            self.n_reads = global_id
            self.metadata_format = metadata_writer.finalize()

        for remap in self.study_remaps.values():
            if hasattr(remap, 'flush'):
                remap.flush()

        for study in study_inputs:
            observed = self.study_remap_counts[study.study_id]
            if observed != study.n_reads:
                raise RuntimeError(
                    f"{study.study_id}: sequence file/remap count mismatch; "
                    f"metadata n_reads={study.n_reads}, observed={observed}"
                )

    def build_global_vocabulary(self, study_inputs: List[StudyInput],
                                temp_dir: Path):
        """
        Build global vocabulary using tournament merge.

        Args:
            study_inputs: Study matrix inputs with paths and metadata
            temp_dir: Temporary directory for merge files
        """
        partition_str = f" [{self.partition}]" if self.partition else ""
        print(f"\nPhase 1: Building global vocabulary{partition_str} via k-way stream", file=sys.stderr)

        self._prepare_study_remaps(study_inputs, temp_dir)
        self._stream_global_vocabulary(study_inputs, temp_dir)

        print(f"  Global vocabulary: {self.n_reads:,} unique sequences", file=sys.stderr)
        print(f"  Total samples: {len(self.all_samples):,}", file=sys.stderr)

    def build_global_matrix(self, study_matrices: Dict[str, Path]) -> Path:
        """
        Build the global matrix by merging all study matrices.
        """
        n_reads = self.n_reads
        n_samples = len(self.all_samples)

        partition_str = f" [{self.partition}]" if self.partition else ""
        print(f"\nPhase 2: Building global matrix{partition_str}: {n_reads:,} reads × {n_samples:,} samples",
              file=sys.stderr)

        if HAVE_ZARR:
            return self._build_zarr_matrix(study_matrices, n_reads, n_samples)
        else:
            return self._build_npz_matrix(study_matrices, n_reads, n_samples)

    def _build_zarr_matrix(self, study_matrices: Dict[str, Path],
                           n_reads: int, n_samples: int) -> Path:
        """Build chunked Zarr matrix."""
        matrix_path = self.output_dir / f'{self.prefix}_matrix.zarr'

        # Initialize Zarr array - handle both v2 and v3 APIs
        zarr_version = tuple(int(x) for x in zarr.__version__.split('.')[:2])

        if zarr_version >= (3, 0):
            # Zarr v3 API - create array directly, then open group for metadata
            matrix_path.mkdir(parents=True, exist_ok=True)
            counts_path = matrix_path / 'counts'
            counts = zarr.open_array(
                str(counts_path),
                mode='w',
                shape=(n_reads, n_samples),
                chunks=(self.chunk_size, n_samples),
                dtype='uint32',
                fill_value=0
            )
            # Open the parent group for metadata
            root = zarr.open_group(str(matrix_path), mode='a')
        else:
            # Zarr v2 API
            blosc_compressor = numcodecs.Blosc(cname='zstd', clevel=3, shuffle=2)
            store = zarr.DirectoryStore(str(matrix_path))
            root = zarr.group(store=store, overwrite=True)
            counts = root.create_dataset(
                'counts',
                shape=(n_reads, n_samples),
                chunks=(self.chunk_size, n_samples),
                dtype='uint32',
                compressor=blosc_compressor,
                fill_value=0
            )

        # Store metadata (attrs API works in both v2 and v3)
        root.attrs['n_reads'] = n_reads
        root.attrs['n_samples'] = n_samples
        root.attrs['samples'] = self.all_samples
        root.attrs['chunk_size'] = self.chunk_size
        root.attrs['partition'] = self.partition

        # Process each study
        for study_id, matrix_path_study in study_matrices.items():
            print(f"  Merging study: {study_id}", file=sys.stderr)

            study_matrix = sp.load_npz(matrix_path_study).tocsr()
            remap = self.study_remaps[study_id]
            sample_offset = self.study_sample_offsets[study_id]
            if study_matrix.shape[0] != len(remap):
                raise RuntimeError(
                    f"{study_id}: matrix row count ({study_matrix.shape[0]:,}) does not match "
                    f"sequence metadata/remap count ({len(remap):,})"
                )

            # Study local rows and global rows are both sorted by sequence, so
            # remapped global chunks are non-decreasing and can be flushed
            # incrementally without a full COO row array.
            current_chunk_idx: Optional[int] = None
            current_entries: List[Tuple[int, int, int]] = []

            def flush_entries(chunk_idx: int, entries: List[Tuple[int, int, int]]):
                if not entries:
                    return
                chunk_start = chunk_idx * self.chunk_size
                chunk_end = min(chunk_start + self.chunk_size, n_reads)

                chunk = counts[chunk_start:chunk_end, :]
                for global_row, global_col, val in entries:
                    local_row = global_row - chunk_start
                    chunk[local_row, global_col] = val
                counts[chunk_start:chunk_end, :] = chunk

            for local_row in range(study_matrix.shape[0]):
                global_row = int(remap[local_row])
                chunk_idx = global_row // self.chunk_size

                if current_chunk_idx is None:
                    current_chunk_idx = chunk_idx
                elif chunk_idx != current_chunk_idx:
                    if chunk_idx < current_chunk_idx:
                        raise RuntimeError(
                            f"Global row mapping for {study_id} is not monotonic; "
                            "study sequences may not be sorted"
                        )
                    flush_entries(current_chunk_idx, current_entries)
                    current_entries = []
                    current_chunk_idx = chunk_idx

                row_start = study_matrix.indptr[local_row]
                row_end = study_matrix.indptr[local_row + 1]
                for data_idx in range(row_start, row_end):
                    global_col = int(study_matrix.indices[data_idx]) + sample_offset
                    current_entries.append((
                        global_row,
                        global_col,
                        int(study_matrix.data[data_idx]),
                    ))

            if current_chunk_idx is not None:
                flush_entries(current_chunk_idx, current_entries)

        return matrix_path

    def _build_npz_matrix(self, study_matrices: Dict[str, Path],
                          n_reads: int, n_samples: int) -> Path:
        """Fallback: build sparse NPZ matrix."""
        matrix_path = self.output_dir / f'{self.prefix}_matrix.npz'

        # Collect all data
        all_rows = []
        all_cols = []
        all_data = []

        for study_id, matrix_path_study in study_matrices.items():
            print(f"  Merging study: {study_id}", file=sys.stderr)

            study_matrix = sp.load_npz(matrix_path_study)
            remap = self.study_remaps[study_id]
            sample_offset = self.study_sample_offsets[study_id]

            coo = study_matrix.tocoo()
            missing_remaps = 0

            for i in range(len(coo.data)):
                local_row_int = int(coo.row[i])
                if local_row_int >= len(remap):
                    missing_remaps += 1
                    continue
                global_row = int(remap[local_row_int])
                global_col = coo.col[i] + sample_offset

                all_rows.append(global_row)
                all_cols.append(global_col)
                all_data.append(coo.data[i])

            if missing_remaps:
                raise RuntimeError(
                    f"{study_id}: {missing_remaps:,} matrix rows were missing from global remap"
                )

        # Build global matrix
        global_matrix = sp.coo_matrix(
            (all_data, (all_rows, all_cols)),
            shape=(n_reads, n_samples),
            dtype=np.uint32
        ).tocsr()

        sp.save_npz(matrix_path, global_matrix)

        # Save sample list separately
        with open(self.output_dir / f'{self.prefix}_samples.json', 'w') as f:
            json.dump(self.all_samples, f)

        return matrix_path

    def save_outputs(self) -> Dict[str, Path]:
        """Save all outputs and return paths."""
        n_reads = self.n_reads

        print(f"Unique reads FASTA: {self.fasta_path}", file=sys.stderr)
        print(f"Read metadata: {self.metadata_path}", file=sys.stderr)

        # Write config
        config = {
            'version': '1.0',
            'matrix_format': 'zarr' if HAVE_ZARR else 'npz',
            'metadata_format': self.metadata_format,
            'chunk_size': self.chunk_size,
            'metadata_shard_rows': self.metadata_shard_rows,
            'write_fasta': self.write_fasta,
            'partition': self.partition,
            'n_reads': n_reads,
            'n_samples': len(self.all_samples),
            'n_studies': len(self.studies),
            'studies': self.studies
        }

        config_path = self.output_dir / f'{self.prefix}_config.json'
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)

        with open(self.output_dir / f'{self.prefix}_samples.json', 'w') as f:
            json.dump(self.all_samples, f)

        return {
            'fasta': self.fasta_path,
            'metadata': self.metadata_path,
            'config': config_path
        }

    def close_remaps(self):
        """Flush and release memmap handles before temporary directory cleanup."""
        for remap in self.study_remaps.values():
            if hasattr(remap, 'flush'):
                remap.flush()
        self.study_remaps.clear()


def load_study_metadata(metadata_path: Path, study_id: str) -> Tuple[int, List[str]]:
    """Load the small per-study metadata JSON needed for global matrix layout."""
    with open(metadata_path, 'r') as handle:
        metadata = json.load(handle)

    if 'n_reads' not in metadata:
        raise ValueError(f"{metadata_path} is missing required field 'n_reads'")
    if 'samples' not in metadata:
        raise ValueError(f"{metadata_path} is missing required field 'samples'")

    n_reads = int(metadata['n_reads'])
    samples = list(metadata['samples'])
    if n_reads < 0:
        raise ValueError(f"{metadata_path} has negative n_reads for {study_id}: {n_reads}")
    return n_reads, samples


def strip_suffix(value: str, suffix: str) -> str:
    if not value.endswith(suffix):
        raise ValueError(f"Expected {value!r} to end with {suffix!r}")
    return value[:-len(suffix)]


def make_tempdir(prefix: str, directory: Path):
    """Create a TemporaryDirectory with cleanup-error tolerance when available."""
    try:
        return tempfile.TemporaryDirectory(
            prefix=prefix,
            dir=directory,
            ignore_cleanup_errors=True,
        )
    except TypeError:
        return tempfile.TemporaryDirectory(prefix=prefix, dir=directory)


def main():
    parser = argparse.ArgumentParser(
        description='Merge study matrices into global Zarr matrix using streaming merge'
    )
    parser.add_argument(
        '--study-dirs',
        nargs='+',
        type=Path,
        required=True,
        help='Directories containing study matrix outputs'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for global matrix'
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=10000,
        help='Zarr chunk size (rows per chunk, default: 10000)'
    )
    parser.add_argument(
        '--metadata-shard-rows',
        type=int,
        default=100_000_000,
        help='Maximum rows per metadata parquet shard (default: 100000000)'
    )
    parser.add_argument(
        '--partition',
        type=str,
        default=None,
        help='Partition identifier (e.g., AA, AC, ..., TT). Looks for {study_id}.{partition}_* files.'
    )
    parser.add_argument(
        '--no-write-fasta',
        action='store_false',
        dest='write_fasta',
        help='Create an empty FASTA placeholder instead of writing global read sequences'
    )

    args = parser.parse_args()

    partition_str = f" [{args.partition}]" if args.partition else ""
    print(f"Merging {len(args.study_dirs)} studies{partition_str} using streaming merge", file=sys.stderr)

    # Discover all study files
    study_inputs: List[StudyInput] = []
    study_matrices: Dict[str, Path] = {}

    for study_dir in args.study_dirs:
        # Find study files - look for metadata to determine study_id without
        # loading the large vocab pickle.
        if args.partition:
            metadata_files = list(study_dir.glob(f'*.{args.partition}_metadata.json'))
        else:
            metadata_files = list(study_dir.glob('*_metadata.json'))
            # Filter out partitioned files
            metadata_files = [f for f in metadata_files if not any(
                f'.{di}_metadata.json' in str(f) for di in
                [f"{a}{b}" for a in "ACGT" for b in "ACGT"]
            )]

        if not metadata_files:
            print(f"Warning: No metadata file in {study_dir} for partition={args.partition}, skipping", file=sys.stderr)
            continue

        # Extract study_id from filename
        metadata_path = metadata_files[0]
        stem = metadata_path.stem  # e.g., "PRJNA123.AA_metadata" or "PRJNA123_metadata"
        if args.partition:
            # Remove ".{partition}_metadata" suffix
            study_id = strip_suffix(stem, f'.{args.partition}_metadata')
            prefix = f'{study_id}.{args.partition}'
        else:
            study_id = strip_suffix(stem, '_metadata')
            prefix = study_id

        seq_path = study_dir / f'{prefix}_sequences.txt.gz'
        matrix_path = study_dir / f'{prefix}_matrix.npz'
        metadata_path = study_dir / f'{prefix}_metadata.json'

        if not all(p.exists() for p in [metadata_path, seq_path, matrix_path]):
            print(f"Warning: Incomplete study files for {study_id}, skipping", file=sys.stderr)
            continue

        n_reads, samples = load_study_metadata(metadata_path, study_id)
        study_inputs.append(StudyInput(
            study_id=study_id,
            sequences_path=seq_path,
            matrix_path=matrix_path,
            metadata_path=metadata_path,
            n_reads=n_reads,
            samples=samples,
        ))
        study_matrices[study_id] = matrix_path

    if not study_inputs:
        print("Error: No valid study inputs found", file=sys.stderr)
        sys.exit(1)

    # Create merger and run
    merger = GlobalMatrixMerger(
        args.output_dir,
        chunk_size=args.chunk_size,
        partition=args.partition,
        metadata_shard_rows=args.metadata_shard_rows,
        write_fasta=args.write_fasta,
    )

    with make_tempdir("global_merge_", args.output_dir) as temp_dir:
        try:
            # Phase 1: Build global vocabulary via streaming merge
            merger.build_global_vocabulary(study_inputs, Path(temp_dir))

            # Phase 2: Build global matrix
            matrix_path = merger.build_global_matrix(study_matrices)
        finally:
            merger.close_remaps()

    # Phase 3: Save outputs
    print("\nPhase 3: Saving outputs", file=sys.stderr)
    outputs = merger.save_outputs()

    print(f"\nGlobal matrix: {matrix_path}", file=sys.stderr)
    print(f"Unique reads FASTA: {outputs['fasta']}", file=sys.stderr)
    print(f"Read metadata: {outputs['metadata']}", file=sys.stderr)
    print(f"Config: {outputs['config']}", file=sys.stderr)


if __name__ == '__main__':
    main()
