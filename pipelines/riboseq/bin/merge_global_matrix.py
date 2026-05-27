#!/usr/bin/env python3
"""
Merge study matrices into a global count store using streaming merge.

Takes all study matrix outputs (optionally for a single partition) and produces:
1. Global matrix/count store
2. Global unique reads FASTA (for single-pass alignment)
3. Metadata parquet dataset (read_id, length, etc.)

Uses k-way sequence streaming for O(n log k) complexity.
When used with partitioned inputs, each invocation handles one dinucleotide partition.
The default global store is sparse partitioned Parquet; dense Zarr/NPZ output is
kept as an explicit compatibility mode for small runs.
"""

import argparse
import contextlib
import gzip
import heapq
import json
import os
import shutil
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
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


class SparseParquetFactWriter:
    """Write sparse global counts as partitioned parquet fact shards."""

    def __init__(
        self,
        output_path: Path,
        temp_dir: Path,
        shard_rows: int,
        read_bucket_size: int,
    ):
        if not HAVE_POLARS:
            raise RuntimeError("polars is required for --matrix-format sparse-parquet")
        self.output_path = Path(output_path)
        self.staging_path = Path(temp_dir) / f"{self.output_path.name}.staging"
        self.lock_path = self.output_path.with_suffix(self.output_path.suffix + ".lock")
        self.shard_rows = max(1, shard_rows)
        self.read_bucket_size = max(1, read_bucket_size)
        self.current_handle = None
        self.current_path: Optional[Path] = None
        self.current_read_bucket: Optional[int] = None
        self.rows_in_shard = 0
        self.total_rows = 0
        self.part_index = 0
        self.parts: List[Dict[str, Union[str, int]]] = []

    def __enter__(self):
        self._acquire_lock()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        if exc_type is not None:
            self.abort()
        self._release_lock()

    def _acquire_lock(self):
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as exc:
            raise RuntimeError(
                f"Refusing to write sparse store while lock exists: {self.lock_path}"
            ) from exc
        with os.fdopen(fd, "w") as handle:
            handle.write(f"pid={os.getpid()}\n")

        if self.staging_path.exists():
            shutil.rmtree(self.staging_path)
        self.staging_path.mkdir(parents=True)

    def _release_lock(self):
        try:
            self.lock_path.unlink()
        except FileNotFoundError:
            pass

    def _read_bucket(self, read_id: int) -> int:
        return int(read_id) // self.read_bucket_size

    def _open_next_shard(self, read_bucket: int):
        self._convert_current_shard()
        bucket_dir = self.staging_path / f"read_bucket={read_bucket:06d}"
        bucket_dir.mkdir(parents=True, exist_ok=True)
        self.current_path = bucket_dir / f"part-{self.part_index:05d}.tsv"
        self.current_handle = open(self.current_path, "w")
        self.current_handle.write(
            "read_bucket\tread_id\tsample_id\tstudy_id_int\tcount\n"
        )
        self.current_read_bucket = read_bucket
        self.rows_in_shard = 0
        self.part_index += 1

    def _convert_current_shard(self):
        if self.current_path is None:
            return
        self.close()
        if self.rows_in_shard == 0:
            try:
                self.current_path.unlink()
            except FileNotFoundError:
                pass
            self.current_path = None
            return

        parquet_path = self.current_path.with_suffix(".parquet")
        self._scan_tsv(self.current_path).sink_parquet(
            parquet_path,
            compression="zstd",
        )
        try:
            self.current_path.unlink()
        except FileNotFoundError:
            pass
        self.parts.append({
            "path": str(parquet_path.relative_to(self.staging_path)),
            "rows": self.rows_in_shard,
            "read_bucket": self.current_read_bucket if self.current_read_bucket is not None else -1,
        })
        self.current_path = None
        self.rows_in_shard = 0

    def _scan_tsv(self, tsv_path: Path):
        scan_kwargs = {"separator": "\t"}
        schema = {
            "read_bucket": pl.UInt32,
            "read_id": pl.UInt64,
            "sample_id": pl.UInt32,
            "study_id_int": pl.UInt32,
            "count": pl.UInt32,
        }
        try:
            return pl.scan_csv(tsv_path, **scan_kwargs, schema_overrides=schema)
        except TypeError:
            return pl.scan_csv(tsv_path, **scan_kwargs, dtypes=schema)

    def write(
        self,
        study_id: str,
        read_id: int,
        sample_id: int,
        study_id_int: int,
        count: int,
    ):
        if count < 0 or count > np.iinfo(np.uint32).max:
            raise ValueError(f"Count out of uint32 range for {study_id}: {count}")
        read_bucket = self._read_bucket(read_id)
        if (
            self.current_handle is None
            or self.current_read_bucket != read_bucket
            or self.rows_in_shard >= self.shard_rows
        ):
            self._open_next_shard(read_bucket)

        self.current_handle.write(
            f"{read_bucket}\t{read_id}\t{sample_id}\t{study_id_int}\t{count}\n"
        )
        self.rows_in_shard += 1
        self.total_rows += 1

    def close(self):
        if self.current_handle is not None:
            self.current_handle.close()
            self.current_handle = None

    def commit(self):
        self._convert_current_shard()
        if self.output_path.exists():
            if self.output_path.is_dir():
                shutil.rmtree(self.output_path)
            else:
                self.output_path.unlink()
        shutil.move(str(self.staging_path), self.output_path)

    def abort(self):
        self.close()
        if self.staging_path.exists():
            shutil.rmtree(self.staging_path)


def latest_generation(manifest: Dict) -> Dict:
    generations = manifest.get("generations") or []
    if generations:
        return sorted(generations, key=lambda item: int(item["generation"]))[-1]
    return {
        "generation": int(manifest.get("generation", 1)),
        "path": manifest.get("path", "global_matrix.parquet"),
        "parents": manifest.get("parents", []),
        "nnz": manifest.get("nnz", 0),
        "parts": manifest.get("parts", []),
    }


def load_sparse_manifest(path: Path) -> Tuple[Path, Dict]:
    """Load a sparse matrix manifest from an output dir or manifest path."""
    path = Path(path)
    if path.is_dir():
        manifest_path = path / "global_matrix_manifest.json"
    else:
        manifest_path = path
    with open(manifest_path, "r") as handle:
        manifest = json.load(handle)
    if manifest.get("matrix_format") != "sparse-parquet":
        raise ValueError(f"{manifest_path} is not a sparse-parquet manifest")
    return manifest_path, manifest


class BufferedNpyReader:
    """Stream primitive values from an unpickled .npy member."""

    def __init__(self, handle, buffer_values: int = 1_000_000):
        version = np.lib.format.read_magic(handle)
        shape, fortran_order, dtype = np.lib.format._read_array_header(
            handle,
            version,
            max_header_size=10000,
        )
        if fortran_order:
            raise ValueError("Fortran-order .npy arrays are not supported for CSR streaming")
        if dtype.hasobject:
            raise ValueError("Object .npy arrays are not supported for CSR streaming")

        self.handle = handle
        self.shape = tuple(shape)
        self.dtype = np.dtype(dtype)
        self.size = int(np.prod(self.shape, dtype=np.int64))
        self.buffer_values = max(1, buffer_values)
        self.buffer = np.asarray([], dtype=self.dtype)
        self.buffer_pos = 0
        self.values_read = 0

    def _fill(self):
        remaining = self.size - self.values_read
        if remaining <= 0:
            self.buffer = np.asarray([], dtype=self.dtype)
            self.buffer_pos = 0
            return
        n_values = min(self.buffer_values, remaining)
        raw = self.handle.read(n_values * self.dtype.itemsize)
        expected = n_values * self.dtype.itemsize
        if len(raw) != expected:
            raise EOFError(f"Unexpected EOF while reading .npy payload: expected {expected}, got {len(raw)}")
        self.buffer = np.frombuffer(raw, dtype=self.dtype)
        self.buffer_pos = 0
        self.values_read += n_values

    def read_values(self, n_values: int):
        if n_values < 0:
            raise ValueError(f"Cannot read negative number of values: {n_values}")
        out = np.empty(n_values, dtype=self.dtype)
        out_pos = 0
        while out_pos < n_values:
            if self.buffer_pos >= len(self.buffer):
                self._fill()
                if len(self.buffer) == 0:
                    raise EOFError(
                        f"Unexpected EOF after {out_pos} of {n_values} requested values"
                    )
            available = len(self.buffer) - self.buffer_pos
            take = min(available, n_values - out_pos)
            out[out_pos:out_pos + take] = self.buffer[self.buffer_pos:self.buffer_pos + take]
            self.buffer_pos += take
            out_pos += take
        return out


def read_npz_scalar_array(npz_path: Path, name: str):
    with np.load(npz_path) as loaded:
        if name not in loaded:
            raise ValueError(f"{npz_path} is missing required array {name!r}")
        return loaded[name]


def iter_csr_npz_rows(npz_path: Path, buffer_values: int = 1_000_000):
    """
    Stream rows from a SciPy-compatible CSR .npz without materializing it.

    Yields:
        local_row, sample_indices, counts
    """
    shape = tuple(int(value) for value in read_npz_scalar_array(npz_path, "shape"))
    if len(shape) != 2:
        raise ValueError(f"{npz_path} has invalid CSR shape: {shape}")

    with zipfile.ZipFile(npz_path, "r") as zf:
        required = {"data.npy", "indices.npy", "indptr.npy"}
        missing = required.difference(zf.namelist())
        if missing:
            raise ValueError(f"{npz_path} is missing CSR arrays: {sorted(missing)}")

        with zf.open("indptr.npy") as indptr_handle, \
                zf.open("indices.npy") as indices_handle, \
                zf.open("data.npy") as data_handle:
            indptr = BufferedNpyReader(indptr_handle, buffer_values=buffer_values)
            indices = BufferedNpyReader(indices_handle, buffer_values=buffer_values)
            data = BufferedNpyReader(data_handle, buffer_values=buffer_values)

            if indptr.size != shape[0] + 1:
                raise ValueError(
                    f"{npz_path} indptr length {indptr.size:,} does not match rows {shape[0]:,}"
                )
            if indices.size != data.size:
                raise ValueError(
                    f"{npz_path} indices length {indices.size:,} does not match data length {data.size:,}"
                )

            previous = int(indptr.read_values(1)[0])
            if previous != 0:
                raise ValueError(f"{npz_path} CSR indptr starts at {previous}, expected 0")

            for local_row in range(shape[0]):
                current = int(indptr.read_values(1)[0])
                if current < previous:
                    raise ValueError(f"{npz_path} CSR indptr decreases at row {local_row}")
                row_nnz = current - previous
                yield (
                    local_row,
                    indices.read_values(row_nnz),
                    data.read_values(row_nnz),
                )
                previous = current

            if previous != data.size:
                raise ValueError(
                    f"{npz_path} CSR indptr ends at {previous:,}, expected nnz {data.size:,}"
                )


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
        sparse_shard_rows: int = 5_000_000,
        sparse_read_bucket_size: int = 100_000,
        write_fasta: bool = True,
        matrix_format: str = "sparse-parquet",
        allow_dense: bool = False,
        dense_chunk_byte_limit: int = 128 * 1024 * 1024,
        append_to: Optional[Path] = None,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.chunk_size = chunk_size
        self.partition = partition  # e.g., "AA", "AC", etc. or None
        self.metadata_shard_rows = metadata_shard_rows
        self.sparse_shard_rows = sparse_shard_rows
        self.sparse_read_bucket_size = sparse_read_bucket_size
        self.metadata_format = 'file'
        self.write_fasta = write_fasta
        self.matrix_format = matrix_format
        self.allow_dense = allow_dense
        self.dense_chunk_byte_limit = dense_chunk_byte_limit
        self.matrix_backend = matrix_format
        self.matrix_manifest: Optional[Dict[str, Union[str, int, List, Dict]]] = None
        self.append_to = Path(append_to) if append_to else None
        self.previous_manifest_path: Optional[Path] = None
        self.previous_manifest: Optional[Dict] = None
        self.generation = 1
        self.parent_generations: List[int] = []
        self.existing_reads: Dict[str, int] = {}
        self.existing_read_rows = []
        self.new_read_rows = []
        self.base_n_reads = 0
        self.base_samples = []
        self.base_studies = []
        self.base_nnz = 0

        # Output file prefix
        self.prefix = f"global.{partition}" if partition else "global"

        # Per-study info
        self.study_remaps: Dict[str, np.ndarray] = {}  # study -> disk-backed local_id -> global_id
        self.study_n_reads: Dict[str, int] = {}
        self.study_remap_counts: Dict[str, int] = {}
        self.all_samples: List[str] = []
        self.study_sample_offsets: Dict[str, int] = {}  # study -> sample column offset
        self.study_id_ints: Dict[str, int] = {}
        self.studies: List[str] = []
        self.n_reads = 0

        self.fasta_path = self.output_dir / f'{self.prefix}_reads.fasta'
        self.metadata_path = self.output_dir / f'{self.prefix}_metadata.parquet'

        if self.append_to:
            self._load_previous_sparse_store(self.append_to)

    def _load_previous_sparse_store(self, previous: Path):
        """Load stable ID state from an existing sparse global store."""
        if self.matrix_format != "sparse-parquet":
            raise ValueError("--append-to is only supported for --matrix-format sparse-parquet")
        if not HAVE_POLARS:
            raise RuntimeError("polars is required for sparse append mode")

        manifest_path, manifest = load_sparse_manifest(previous)
        self.previous_manifest_path = manifest_path
        self.previous_manifest = manifest

        latest = latest_generation(manifest)
        self.generation = int(latest["generation"]) + 1
        self.parent_generations = [int(latest["generation"])]

        reads_path = manifest_path.parent / manifest.get("global_reads", "global_reads.parquet")
        if not reads_path.exists():
            raise FileNotFoundError(
                f"Append mode requires canonical read lookup table: {reads_path}"
            )
        reads = pl.read_parquet(reads_path).select("read_id", "sequence", "length")
        self.existing_read_rows = reads.sort("read_id").to_dicts()
        self.existing_reads = {
            row["sequence"]: int(row["read_id"])
            for row in self.existing_read_rows
        }
        self.base_n_reads = len(self.existing_reads)

        samples_path = manifest_path.parent / manifest.get("lookup_tables", {}).get("samples", "global_samples.parquet")
        studies_path = manifest_path.parent / manifest.get("lookup_tables", {}).get("studies", "global_studies.parquet")
        if samples_path.exists():
            self.base_samples = pl.read_parquet(samples_path).sort("sample_id").to_dicts()
            self.all_samples = [row["sample_name"] for row in self.base_samples]
        if studies_path.exists():
            self.base_studies = pl.read_parquet(studies_path).sort("study_id_int").to_dicts()
            self.studies = [row["study_id"] for row in self.base_studies]
            self.study_id_ints = {
                row["study_id"]: int(row["study_id_int"])
                for row in self.base_studies
            }
            self.study_n_reads = {
                row["study_id"]: int(row["n_reads"])
                for row in self.base_studies
            }

        self.base_nnz = sum(int(gen.get("nnz", 0)) for gen in manifest.get("generations", [latest]))

    def _prepare_study_remaps(self, study_inputs: List[StudyInput], temp_dir: Path):
        """Create disk-backed remap arrays and record sample offsets."""
        base_study_count = len(self.studies)
        for study_offset, study in enumerate(study_inputs):
            study_idx = base_study_count + study_offset
            if study.study_id in self.study_id_ints:
                raise ValueError(
                    f"Study ID already exists in this global store: {study.study_id}"
                )
            duplicate_samples = set(study.samples).intersection(self.all_samples)
            if duplicate_samples:
                examples = ', '.join(sorted(duplicate_samples)[:5])
                raise ValueError(
                    f"Duplicate sample IDs across study matrices before {study.study_id}: {examples}"
                )

            self.study_sample_offsets[study.study_id] = len(self.all_samples)
            self.study_id_ints[study.study_id] = study_idx
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

            next_global_id = self.base_n_reads
            while heap:
                seq = heap[0][0]
                occurrences: List[Tuple[int, int]] = []

                while heap and heap[0][0] == seq:
                    _, study_index, local_id = heapq.heappop(heap)
                    occurrences.append((study_index, local_id))
                    advance(study_index)

                if seq in self.existing_reads:
                    global_id = self.existing_reads[seq]
                else:
                    global_id = next_global_id
                    next_global_id += 1
                    self.existing_reads[seq] = global_id
                    self.new_read_rows.append({
                        "read_id": global_id,
                        "sequence": seq,
                        "length": len(seq),
                    })

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

            self.n_reads = next_global_id
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

        if self.matrix_format == "sparse-parquet":
            return self._build_sparse_parquet_matrix(study_matrices, n_reads, n_samples)

        self._guard_dense_matrix(n_reads, n_samples)
        if HAVE_ZARR:
            self.matrix_backend = "zarr"
            return self._build_zarr_matrix(study_matrices, n_reads, n_samples)
        else:
            self.matrix_backend = "npz"
            return self._build_npz_matrix(study_matrices, n_reads, n_samples)

    def _guard_dense_matrix(self, n_reads: int, n_samples: int):
        """Fail fast for dense chunk shapes that are too large for production-scale runs."""
        if self.allow_dense:
            return
        estimated_chunk_bytes = min(self.chunk_size, max(n_reads, 1)) * max(n_samples, 1) * np.dtype(np.uint32).itemsize
        if estimated_chunk_bytes > self.dense_chunk_byte_limit:
            raise RuntimeError(
                "Dense global matrix output is unsafe for this run: "
                f"estimated row chunk is {estimated_chunk_bytes:,} bytes "
                f"({min(self.chunk_size, max(n_reads, 1)):,} rows × {n_samples:,} samples × uint32). "
                "Use --matrix-format sparse-parquet for the canonical sparse store, "
                "or pass --allow-dense only for a deliberately small/controlled run."
            )

    def _write_sparse_manifest(
        self,
        generation_path: Path,
        writer: SparseParquetFactWriter,
        n_reads: int,
        n_samples: int,
    ):
        generation_entry = {
            "generation": self.generation,
            "path": str(generation_path.relative_to(self.output_dir)),
            "parents": self.parent_generations,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "nnz": writer.total_rows,
            "parts": writer.parts,
            "n_new_reads": len(self.new_read_rows),
            "n_new_samples": n_samples - len(self.base_samples),
            "n_new_studies": len(self.studies) - len(self.base_studies),
        }
        previous_generations = []
        if self.previous_manifest:
            previous_generations = list(self.previous_manifest.get("generations", []))
            if not previous_generations:
                previous_latest = latest_generation(self.previous_manifest)
                previous_generations = [previous_latest]
        generations = previous_generations + [generation_entry]
        total_nnz = sum(int(item.get("nnz", 0)) for item in generations)

        manifest = {
            "version": "1.0",
            "matrix_format": "sparse-parquet",
            "status": "committed",
            "generation": self.generation,
            "parents": self.parent_generations,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "path": generation_entry["path"],
            "global_counts": "global_counts",
            "global_reads": "global_reads.parquet",
            "tombstones_path": "global_tombstones.parquet",
            "retained_counts": "global_retained_counts.parquet",
            "partition": self.partition,
            "schema": {
                "read_bucket": "uint32",
                "read_id": "uint64",
                "sample_id": "uint32",
                "study_id_int": "uint32",
                "count": "uint32",
            },
            "lookup_tables": {
                "samples": f"{self.prefix}_samples.parquet",
                "studies": f"{self.prefix}_studies.parquet",
            },
            "partitioning": ["generation", "read_bucket"],
            "read_bucket_size": self.sparse_read_bucket_size,
            "n_reads": n_reads,
            "n_samples": n_samples,
            "n_studies": len(self.studies),
            "nnz": total_nnz,
            "parts": writer.parts,
            "generations": generations,
            "append_protocol": {
                "mode": "generation",
                "description": "Append by writing a new sparse generation with this manifest in parents; do not destructively filter or rewrite committed facts.",
            },
            "commit_protocol": {
                "lock": f"{generation_path.name}.lock",
                "staging": f"{generation_path.name}.staging",
                "commit": "write parts in staging, atomically move staging into place, write manifest, release lock",
            },
        }
        manifest_path = self.output_dir / f"{self.prefix}_matrix_manifest.json"
        with open(manifest_path, "w") as handle:
            json.dump(manifest, handle, indent=2)
        self.matrix_manifest = manifest

    def _build_sparse_parquet_matrix(
        self,
        study_matrices: Dict[str, Path],
        n_reads: int,
        n_samples: int,
    ) -> Path:
        """Build partitioned parquet sparse fact table without dense row chunks."""
        counts_root = self.output_dir / f"{self.prefix}_counts"
        matrix_path = counts_root / f"generation={self.generation:06d}"
        compat_path = self.output_dir / f"{self.prefix}_matrix.parquet"
        print("  Writing sparse parquet facts", file=sys.stderr)

        if self.previous_manifest_path and self.previous_manifest_path.parent != self.output_dir:
            for generation in self.previous_manifest.get("generations", [latest_generation(self.previous_manifest)]):
                source = self.previous_manifest_path.parent / generation["path"]
                destination = self.output_dir / generation["path"]
                if destination.exists():
                    continue
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copytree(source, destination)

        with make_tempdir("global_sparse_", self.output_dir) as temp_dir:
            with SparseParquetFactWriter(
                matrix_path,
                Path(temp_dir),
                self.sparse_shard_rows,
                self.sparse_read_bucket_size,
            ) as writer:
                for study_id, matrix_path_study in study_matrices.items():
                    print(f"  Merging study: {study_id}", file=sys.stderr)

                    remap = self.study_remaps[study_id]
                    sample_offset = self.study_sample_offsets[study_id]
                    shape = tuple(int(value) for value in read_npz_scalar_array(matrix_path_study, "shape"))
                    if shape[0] != len(remap):
                        raise RuntimeError(
                            f"{study_id}: matrix row count ({shape[0]:,}) does not match "
                            f"sequence metadata/remap count ({len(remap):,})"
                        )

                    for local_row, sample_indices, counts in iter_csr_npz_rows(matrix_path_study):
                        global_row = int(remap[local_row])
                        for sample_idx, count in zip(sample_indices, counts):
                            sample_index = int(sample_idx) + sample_offset
                            writer.write(
                                study_id=study_id,
                                read_id=global_row,
                                sample_id=sample_index,
                                study_id_int=self.study_id_ints[study_id],
                                count=int(count),
                            )

                writer.commit()
                self._write_sparse_manifest(matrix_path, writer, n_reads, n_samples)

        if compat_path.exists():
            if compat_path.is_dir():
                shutil.rmtree(compat_path)
            else:
                compat_path.unlink()
        if matrix_path.exists():
            shutil.copytree(matrix_path, compat_path)

        return matrix_path

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

    def _write_sparse_lookup_tables(self):
        """Write compact ID lookup tables for sparse fact rows."""
        if self.matrix_format != "sparse-parquet":
            return
        if not HAVE_POLARS:
            raise RuntimeError("polars is required for sparse lookup tables")

        sample_rows = list(self.base_samples)
        new_studies = [study_id for study_id in self.studies if study_id in self.study_sample_offsets]
        for study_id in new_studies:
            study_id_int = self.study_id_ints[study_id]
            sample_offset = self.study_sample_offsets[study_id]
            ordered_new_studies = sorted(new_studies, key=lambda item: self.study_id_ints[item])
            current_pos = ordered_new_studies.index(study_id)
            if current_pos + 1 < len(ordered_new_studies):
                next_offset = self.study_sample_offsets[ordered_new_studies[current_pos + 1]]
            else:
                next_offset = len(self.all_samples)

            for sample_id in range(sample_offset, next_offset):
                sample_rows.append({
                    "sample_id": sample_id,
                    "sample_name": self.all_samples[sample_id],
                    "study_id_int": study_id_int,
                    "study_id": study_id,
                    "study_sample_index": sample_id - sample_offset,
                })

        pl.DataFrame(
            sample_rows,
            schema={
                "sample_id": pl.UInt32,
                "sample_name": pl.Utf8,
                "study_id_int": pl.UInt32,
                "study_id": pl.Utf8,
                "study_sample_index": pl.UInt32,
            },
        ).write_parquet(self.output_dir / f"{self.prefix}_samples.parquet")

        study_rows = list(self.base_studies)
        study_rows.extend(
            [
                {
                    "study_id_int": self.study_id_ints[study_id],
                    "study_id": study_id,
                    "n_reads": self.study_n_reads[study_id],
                    "sample_offset": self.study_sample_offsets[study_id],
                }
                for study_id in new_studies
            ]
        )
        pl.DataFrame(
            study_rows,
            schema={
                "study_id_int": pl.UInt32,
                "study_id": pl.Utf8,
                "n_reads": pl.UInt64,
                "sample_offset": pl.UInt32,
            },
        ).write_parquet(self.output_dir / f"{self.prefix}_studies.parquet")

    def _write_global_read_tables(self):
        """Write canonical read lookup table and compatibility metadata."""
        if self.matrix_format != "sparse-parquet":
            return
        if not HAVE_POLARS:
            raise RuntimeError("polars is required for sparse read lookup tables")

        rows = list(self.existing_read_rows)
        existing_ids = {int(row["read_id"]) for row in rows}
        rows.extend(row for row in self.new_read_rows if int(row["read_id"]) not in existing_ids)
        rows = sorted(rows, key=lambda row: int(row["read_id"]))

        frame = pl.DataFrame(
            rows,
            schema={
                "read_id": pl.UInt64,
                "sequence": pl.Utf8,
                "length": pl.UInt32,
            },
        ).with_columns(
            (pl.col("read_id") // self.sparse_read_bucket_size)
            .cast(pl.UInt32)
            .alias("read_bucket")
        ).select(
            "read_id",
            "read_bucket",
            "sequence",
            "length",
        )
        frame.write_parquet(self.output_dir / f"{self.prefix}_reads.parquet")
        if self.prefix == "global":
            frame.write_parquet(self.output_dir / "global_reads.parquet")

        if self.metadata_path.exists():
            if self.metadata_path.is_dir():
                shutil.rmtree(self.metadata_path)
            else:
                self.metadata_path.unlink()
        frame.select("read_id", "length").write_parquet(self.metadata_path)
        self.metadata_format = "file"

        if self.write_fasta:
            with open(self.fasta_path, "w") as handle:
                for row in rows:
                    handle.write(f">read_{int(row['read_id'])}\n{row['sequence']}\n")

    def _ensure_tombstones_table(self):
        if self.matrix_format != "sparse-parquet":
            return
        tombstones_path = self.output_dir / f"{self.prefix}_tombstones.parquet"
        if tombstones_path.exists():
            return
        if self.previous_manifest_path:
            previous_tombstones = (
                self.previous_manifest_path.parent
                / self.previous_manifest.get("tombstones_path", "global_tombstones.parquet")
            )
            if previous_tombstones.exists():
                shutil.copy2(previous_tombstones, tombstones_path)
                if self.prefix == "global":
                    shutil.copy2(previous_tombstones, self.output_dir / "global_tombstones.parquet")
                return

        empty = pl.DataFrame(
            {
                "tombstone_id": pl.Series([], dtype=pl.UInt64),
                "active": pl.Series([], dtype=pl.Boolean),
                "created_at": pl.Series([], dtype=pl.Utf8),
                "reason": pl.Series([], dtype=pl.Utf8),
                "read_id": pl.Series([], dtype=pl.UInt64),
                "sample_id": pl.Series([], dtype=pl.Utf8),
                "study_id": pl.Series([], dtype=pl.Utf8),
            }
        )
        empty.write_parquet(tombstones_path)
        if self.prefix == "global":
            empty.write_parquet(self.output_dir / "global_tombstones.parquet")

    def _write_retained_counts_view(self):
        """Materialize a query-friendly retained-count fact view."""
        if self.matrix_format != "sparse-parquet" or not self.matrix_manifest:
            return
        retained_path = self.output_dir / f"{self.prefix}_retained_counts.parquet"
        if retained_path.exists():
            if retained_path.is_dir():
                shutil.rmtree(retained_path)
            else:
                retained_path.unlink()

        paths = [
            str(self.output_dir / generation["path"] / "read_bucket=*" / "*.parquet")
            for generation in self.matrix_manifest.get("generations", [])
        ]
        if not paths:
            return

        scan = pl.concat([pl.scan_parquet(path) for path in paths])
        samples = pl.scan_parquet(str(self.output_dir / f"{self.prefix}_samples.parquet")).select(
            "sample_id",
            "sample_name",
            "study_id",
        )
        scan = scan.join(samples, on="sample_id", how="left")

        tombstones_path = self.output_dir / f"{self.prefix}_tombstones.parquet"
        if tombstones_path.exists():
            tombstones = pl.read_parquet(tombstones_path).filter(pl.col("active") == True).to_dicts()
            expr = None
            for tombstone in tombstones:
                current = None
                if tombstone.get("read_id") is not None:
                    current = pl.col("read_id") == int(tombstone["read_id"])
                if tombstone.get("sample_id") is not None:
                    sample_expr = pl.col("sample_name") == str(tombstone["sample_id"])
                    current = sample_expr if current is None else current & sample_expr
                if tombstone.get("study_id") is not None:
                    study_expr = pl.col("study_id") == str(tombstone["study_id"])
                    current = study_expr if current is None else current & study_expr
                if current is not None:
                    expr = current if expr is None else expr | current
            if expr is not None:
                scan = scan.filter(~expr)

        scan.select("read_bucket", "read_id", "sample_id", "study_id_int", "count").sink_parquet(
            retained_path,
            compression="zstd",
        )
        if self.prefix == "global":
            compat = self.output_dir / "global_retained_counts.parquet"
            if compat != retained_path:
                shutil.copy2(retained_path, compat)

    def save_outputs(self) -> Dict[str, Path]:
        """Save all outputs and return paths."""
        n_reads = self.n_reads

        print(f"Unique reads FASTA: {self.fasta_path}", file=sys.stderr)
        print(f"Read metadata: {self.metadata_path}", file=sys.stderr)
        self._write_global_read_tables()
        self._write_sparse_lookup_tables()
        self._ensure_tombstones_table()
        self._write_retained_counts_view()

        # Write config
        config = {
            'version': '1.0',
            'matrix_format': self.matrix_backend,
            'metadata_format': self.metadata_format,
            'chunk_size': self.chunk_size,
            'metadata_shard_rows': self.metadata_shard_rows,
            'sparse_shard_rows': self.sparse_shard_rows,
            'sparse_read_bucket_size': self.sparse_read_bucket_size,
            'write_fasta': self.write_fasta,
            'partition': self.partition,
            'n_reads': n_reads,
            'n_samples': len(self.all_samples),
            'n_studies': len(self.studies),
            'studies': self.studies,
            'generation': self.generation if self.matrix_manifest else None,
            'global_reads': f'{self.prefix}_reads.parquet' if self.matrix_manifest else None,
            'global_counts': f'{self.prefix}_counts' if self.matrix_manifest else None,
            'tombstones': f'{self.prefix}_tombstones.parquet' if self.matrix_manifest else None,
            'retained_counts': f'{self.prefix}_retained_counts.parquet' if self.matrix_manifest else None,
            'matrix_manifest': f'{self.prefix}_matrix_manifest.json' if self.matrix_manifest else None,
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
        description='Merge study matrices into a global count store using streaming merge'
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
        help='Dense/Zarr chunk size in rows (default: 10000)'
    )
    parser.add_argument(
        '--metadata-shard-rows',
        type=int,
        default=100_000_000,
        help='Maximum rows per metadata parquet shard (default: 100000000)'
    )
    parser.add_argument(
        '--matrix-format',
        choices=['sparse-parquet', 'dense'],
        default='sparse-parquet',
        help='Global matrix output format. sparse-parquet is the scalable default; dense keeps legacy Zarr/NPZ behavior.'
    )
    parser.add_argument(
        '--sparse-shard-rows',
        type=int,
        default=5_000_000,
        help='Maximum nonzero count rows per sparse parquet shard (default: 5000000)'
    )
    parser.add_argument(
        '--sparse-read-bucket-size',
        type=int,
        default=100_000,
        help='Number of global read IDs per sparse parquet read_bucket partition (default: 100000)'
    )
    parser.add_argument(
        '--allow-dense',
        action='store_true',
        help='Allow dense global output even when the estimated dense row chunk exceeds the safety limit'
    )
    parser.add_argument(
        '--dense-chunk-byte-limit',
        type=int,
        default=128 * 1024 * 1024,
        help='Maximum estimated dense row chunk bytes before --matrix-format dense fails without --allow-dense'
    )
    parser.add_argument(
        '--append-to',
        type=Path,
        help='Existing sparse global output directory or manifest to append as a new generation'
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
        sparse_shard_rows=args.sparse_shard_rows,
        sparse_read_bucket_size=args.sparse_read_bucket_size,
        write_fasta=args.write_fasta,
        matrix_format=args.matrix_format,
        allow_dense=args.allow_dense,
        dense_chunk_byte_limit=args.dense_chunk_byte_limit,
        append_to=args.append_to,
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
