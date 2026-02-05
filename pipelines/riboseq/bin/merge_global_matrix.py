#!/usr/bin/env python3
"""
Merge study matrices into a global matrix with Zarr output using tournament-style merge.

Takes all study matrix outputs (optionally for a single partition) and produces:
1. Global Zarr matrix (chunked by reads for efficient locus queries)
2. Global unique reads FASTA (for single-pass alignment)
3. Metadata parquet (read_id, length, etc.)

Uses tournament-style pairwise merging for O(n log k) complexity.
When used with partitioned inputs, each invocation handles one dinucleotide partition.
Designed for 100-250 studies, ~100-200M total unique reads.
"""

import argparse
import gzip
import json
import os
import pickle
import struct
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

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


# =============================================================================
# Streaming merge helpers
# =============================================================================

def read_study_sequences(sequences_path: Path, study_id: str) -> Iterator[Tuple[str, str]]:
    """
    Stream sorted sequences from a study's gzipped sequences file.
    Yields: (sequence, study_id)
    """
    with gzip.open(sequences_path, 'rt') as f:
        for line in f:
            seq = line.strip()
            if seq:
                yield (seq, study_id)


def read_intermediate_file(path: Path) -> Iterator[Tuple[str, List[str]]]:
    """
    Read intermediate merge file.
    Format: sequence\tstudy1,study2,...

    Yields: (sequence, [study_ids])
    """
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            seq = parts[0]
            studies = parts[1].split(',') if len(parts) > 1 and parts[1] else []
            yield seq, studies


def write_intermediate_file(path: Path, data: Iterator[Tuple[str, List[str]]]):
    """
    Write intermediate merge file.
    Format: sequence\tstudy1,study2,...
    """
    with open(path, 'w') as f:
        for seq, studies in data:
            studies_str = ','.join(studies)
            f.write(f"{seq}\t{studies_str}\n")


def merge_two_sources(
    source1: Iterator[Tuple[str, List[str]]],
    source2: Iterator[Tuple[str, List[str]]]
) -> Iterator[Tuple[str, List[str]]]:
    """
    Two-pointer merge of two sorted sources.
    Each source yields (sequence, [study_ids])

    This is the core merge operation - O(n) for two sorted lists.
    """
    item1: Optional[Tuple[str, List[str]]] = None
    item2: Optional[Tuple[str, List[str]]] = None

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
            seq1, studies1 = item1
            seq2, studies2 = item2

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
                # Same sequence - merge study lists
                merged_studies = studies1 + studies2
                yield (seq1, merged_studies)
                try:
                    item1 = next(source1)
                except StopIteration:
                    item1 = None
                try:
                    item2 = next(source2)
                except StopIteration:
                    item2 = None


def study_to_intermediate(sequences_path: Path, study_id: str) -> Iterator[Tuple[str, List[str]]]:
    """Convert a study's sequences to intermediate format."""
    for seq, sid in read_study_sequences(sequences_path, study_id):
        yield (seq, [sid])


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

    def __init__(self, output_dir: Path, chunk_size: int = 10000, partition: str = None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.chunk_size = chunk_size
        self.partition = partition  # e.g., "AA", "AC", etc. or None

        # Output file prefix
        self.prefix = f"global.{partition}" if partition else "global"

        # Global sequence tracking
        self.sequences: List[str] = []  # In sorted order from tournament merge
        self.seq_to_global_id: Dict[str, int] = {}  # sequence -> global_id

        # Per-study info
        self.study_local_to_global: Dict[str, Dict[int, int]] = {}  # study -> {local_id: global_id}
        self.all_samples: List[str] = []
        self.study_sample_offsets: Dict[str, int] = {}  # study -> sample column offset
        self.studies: List[str] = []

        # Study vocabs (loaded once, used for remapping)
        self.study_vocabs: Dict[str, Dict[int, int]] = {}  # study_id -> {seq_hash: local_id}

    def _load_study_vocabs(self, study_inputs: List[Tuple[str, Path, Path, Path]]):
        """Load all study vocabularies and record sample offsets."""
        for study_id, vocab_path, _, _ in study_inputs:
            with open(vocab_path, 'rb') as f:
                vocab_data = pickle.load(f)

            self.study_vocabs[study_id] = vocab_data['seq_to_id']
            study_samples = vocab_data['samples']

            # Record sample offset
            self.study_sample_offsets[study_id] = len(self.all_samples)
            self.all_samples.extend(study_samples)
            self.studies.append(study_id)

            # Initialize remap dict
            self.study_local_to_global[study_id] = {}

    def _process_merged_file(self, merged_file: Path):
        """Process merged file and build global structures."""
        for seq, studies_containing in read_intermediate_file(merged_file):
            global_id = len(self.sequences)
            self.sequences.append(seq)
            self.seq_to_global_id[seq] = global_id

            # Build remapping for each study that has this sequence
            seq_hash = hash_sequence(seq)
            for study_id in studies_containing:
                if study_id in self.study_vocabs:
                    local_id = self.study_vocabs[study_id].get(seq_hash)
                    if local_id is not None:
                        self.study_local_to_global[study_id][local_id] = global_id

    def build_global_vocabulary(self, study_inputs: List[Tuple[str, Path, Path, Path]],
                                temp_dir: Path):
        """
        Build global vocabulary using tournament merge.

        Args:
            study_inputs: List of (study_id, vocab_path, sequences_path, matrix_path)
            temp_dir: Temporary directory for merge files
        """
        partition_str = f" [{self.partition}]" if self.partition else ""
        print(f"\nPhase 1: Building global vocabulary{partition_str} via tournament merge", file=sys.stderr)

        # Load all study vocabs first
        self._load_study_vocabs(study_inputs)

        # Prepare inputs for tournament
        merge_inputs = [(study_id, seq_path) for study_id, _, seq_path, _ in study_inputs]

        merger = TournamentMerger(temp_dir)
        final_file = merger.merge_studies(merge_inputs)

        print(f"\n  Building global ID mappings from merged sequences", file=sys.stderr)
        self._process_merged_file(final_file)
        merger.cleanup(keep_final=False)

        print(f"  Global vocabulary: {len(self.sequences):,} unique sequences", file=sys.stderr)
        print(f"  Total samples: {len(self.all_samples):,}", file=sys.stderr)

    def build_global_matrix(self, study_matrices: Dict[str, Path]) -> Path:
        """
        Build the global matrix by merging all study matrices.
        """
        n_reads = len(self.sequences)
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

            study_matrix = sp.load_npz(matrix_path_study)
            remap = self.study_local_to_global[study_id]
            sample_offset = self.study_sample_offsets[study_id]

            # Convert to COO for iteration
            coo = study_matrix.tocoo()

            # Group by chunk for efficient writing
            chunk_data = {}  # chunk_idx -> list of (global_row, global_col, val)

            for i in range(len(coo.data)):
                local_row = coo.row[i]
                col = coo.col[i]
                val = coo.data[i]

                global_row = remap.get(local_row)
                if global_row is None:
                    continue  # Sequence not in global vocab (shouldn't happen)

                global_col = col + sample_offset
                chunk_idx = global_row // self.chunk_size

                if chunk_idx not in chunk_data:
                    chunk_data[chunk_idx] = []
                chunk_data[chunk_idx].append((global_row, global_col, val))

            # Write chunks
            for chunk_idx, entries in chunk_data.items():
                chunk_start = chunk_idx * self.chunk_size
                chunk_end = min(chunk_start + self.chunk_size, n_reads)

                # Read current chunk
                chunk = counts[chunk_start:chunk_end, :]

                # Apply updates
                for global_row, global_col, val in entries:
                    local_row = global_row - chunk_start
                    chunk[local_row, global_col] = val

                # Write back
                counts[chunk_start:chunk_end, :] = chunk

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
            remap = self.study_local_to_global[study_id]
            sample_offset = self.study_sample_offsets[study_id]

            coo = study_matrix.tocoo()

            for i in range(len(coo.data)):
                global_row = remap.get(coo.row[i])
                if global_row is None:
                    continue
                global_col = coo.col[i] + sample_offset

                all_rows.append(global_row)
                all_cols.append(global_col)
                all_data.append(coo.data[i])

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
        n_reads = len(self.sequences)

        # Write FASTA
        fasta_path = self.output_dir / f'{self.prefix}_reads.fasta'
        print(f"Writing FASTA: {fasta_path}", file=sys.stderr)
        with open(fasta_path, 'w') as fasta:
            for global_id, seq in enumerate(self.sequences):
                fasta.write(f">read_{global_id}\n{seq}\n")

        # Write metadata
        metadata_path = self.output_dir / f'{self.prefix}_metadata.parquet'
        print(f"Writing metadata: {metadata_path}", file=sys.stderr)

        if HAVE_POLARS:
            read_ids = list(range(n_reads))
            lengths = [len(seq) for seq in self.sequences]

            df = pl.DataFrame({
                'read_id': read_ids,
                'length': lengths
            })
            df.write_parquet(metadata_path)
        else:
            # Fallback to JSON
            metadata = [
                {'read_id': i, 'length': len(seq)}
                for i, seq in enumerate(self.sequences)
            ]
            with open(metadata_path.with_suffix('.json'), 'w') as f:
                json.dump(metadata, f)

        # Write config
        config = {
            'version': '1.0',
            'matrix_format': 'zarr' if HAVE_ZARR else 'npz',
            'chunk_size': self.chunk_size,
            'partition': self.partition,
            'n_reads': n_reads,
            'n_samples': len(self.all_samples),
            'n_studies': len(self.studies),
            'studies': self.studies
        }

        config_path = self.output_dir / f'{self.prefix}_config.json'
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)

        return {
            'fasta': fasta_path,
            'metadata': metadata_path,
            'config': config_path
        }


def main():
    parser = argparse.ArgumentParser(
        description='Merge study matrices into global Zarr matrix using tournament merge'
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
        '--partition',
        type=str,
        default=None,
        help='Partition identifier (e.g., AA, AC, ..., TT). Looks for {study_id}.{partition}_* files.'
    )

    args = parser.parse_args()

    partition_str = f" [{args.partition}]" if args.partition else ""
    print(f"Merging {len(args.study_dirs)} studies{partition_str} using tournament merge", file=sys.stderr)

    # Discover all study files
    study_inputs = []
    study_matrices = {}

    # File pattern depends on whether partitioned
    if args.partition:
        suffix = f".{args.partition}"
    else:
        suffix = ""

    for study_dir in args.study_dirs:
        # Find study files - look for vocab to determine study_id
        if args.partition:
            vocab_files = list(study_dir.glob(f'*.{args.partition}_vocab.pkl'))
        else:
            vocab_files = list(study_dir.glob('*_vocab.pkl'))
            # Filter out partitioned files
            vocab_files = [f for f in vocab_files if not any(
                f'.{di}_vocab.pkl' in str(f) for di in
                [f"{a}{b}" for a in "ACGT" for b in "ACGT"]
            )]

        if not vocab_files:
            print(f"Warning: No vocab file in {study_dir} for partition={args.partition}, skipping", file=sys.stderr)
            continue

        # Extract study_id from filename
        vocab_file = vocab_files[0]
        stem = vocab_file.stem  # e.g., "PRJNA123.AA_vocab" or "PRJNA123_vocab"
        if args.partition:
            # Remove ".{partition}_vocab" suffix
            study_id = stem.replace(f'.{args.partition}_vocab', '')
            prefix = f'{study_id}.{args.partition}'
        else:
            study_id = stem.replace('_vocab', '')
            prefix = study_id

        vocab_path = study_dir / f'{prefix}_vocab.pkl'
        seq_path = study_dir / f'{prefix}_sequences.txt.gz'
        matrix_path = study_dir / f'{prefix}_matrix.npz'

        if not all(p.exists() for p in [vocab_path, seq_path, matrix_path]):
            print(f"Warning: Incomplete study files for {study_id}, skipping", file=sys.stderr)
            continue

        study_inputs.append((study_id, vocab_path, seq_path, matrix_path))
        study_matrices[study_id] = matrix_path

    if not study_inputs:
        print("Error: No valid study inputs found", file=sys.stderr)
        sys.exit(1)

    # Create merger and run
    merger = GlobalMatrixMerger(args.output_dir, chunk_size=args.chunk_size, partition=args.partition)

    with tempfile.TemporaryDirectory(prefix="global_merge_") as temp_dir:
        # Phase 1: Build global vocabulary via tournament merge
        merger.build_global_vocabulary(study_inputs, Path(temp_dir))

        # Phase 2: Build global matrix
        matrix_path = merger.build_global_matrix(study_matrices)

    # Phase 3: Save outputs
    print("\nPhase 3: Saving outputs", file=sys.stderr)
    outputs = merger.save_outputs()

    print(f"\nGlobal matrix: {matrix_path}", file=sys.stderr)
    print(f"Unique reads FASTA: {outputs['fasta']}", file=sys.stderr)
    print(f"Read metadata: {outputs['metadata']}", file=sys.stderr)
    print(f"Config: {outputs['config']}", file=sys.stderr)


if __name__ == '__main__':
    main()
