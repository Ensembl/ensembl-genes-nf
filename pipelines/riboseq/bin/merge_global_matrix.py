#!/usr/bin/env python3
"""
Merge study matrices into a global matrix with Zarr output.

Takes all study matrix outputs and produces:
1. Global Zarr matrix (chunked by reads for efficient locus queries)
2. Global unique reads FASTA (for single-pass alignment)
3. Metadata parquet (read_id, length, etc.)

Designed for 100-250 studies, ~100-200M total unique reads.
Uses memory-mapped sequence storage to avoid OOM.
"""

import argparse
import gzip
import json
import mmap
import pickle
import struct
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


class MemoryMappedSequenceStore:
    """
    Memory-efficient sequence storage using memory-mapped files.

    Stores sequences on disk, keeps only hash->id mapping in RAM.
    """

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Memory-mapped file for sequences
        self.seq_file_path = self.output_dir / '_sequences.mmap'
        self.seq_file = open(self.seq_file_path, 'wb+')

        # Index: [(offset, length), ...]
        self.seq_index: List[Tuple[int, int]] = []

        # Hash -> global_id mapping (this stays in RAM)
        self.hash_to_id: Dict[int, int] = {}

        self.next_id = 0

    def add_sequence(self, seq: str, seq_hash: int) -> Tuple[int, bool]:
        """
        Add sequence to store if not already present.

        Returns:
            (global_id, is_new) tuple
        """
        if seq_hash in self.hash_to_id:
            return self.hash_to_id[seq_hash], False

        # Write sequence to mmap file
        offset = self.seq_file.tell()
        seq_bytes = seq.encode('utf-8')

        # Write: 4-byte length + sequence bytes
        self.seq_file.write(struct.pack('I', len(seq_bytes)))
        self.seq_file.write(seq_bytes)

        # Record in index
        self.seq_index.append((offset, len(seq_bytes)))

        # Assign ID
        global_id = self.next_id
        self.hash_to_id[seq_hash] = global_id
        self.next_id += 1

        return global_id, True

    def get_sequence(self, global_id: int) -> str:
        """Retrieve sequence by global ID."""
        offset, length = self.seq_index[global_id]
        self.seq_file.seek(offset)
        _ = struct.unpack('I', self.seq_file.read(4))[0]
        return self.seq_file.read(length).decode('utf-8')

    def finalize(self):
        """Flush and prepare for reading."""
        self.seq_file.flush()

    def write_fasta(self, fasta_path: Path):
        """Write all sequences to FASTA file."""
        self.seq_file.seek(0)

        with open(fasta_path, 'w') as fasta:
            for global_id in range(self.next_id):
                seq = self.get_sequence(global_id)
                fasta.write(f">read_{global_id}\n{seq}\n")

    def write_metadata(self, metadata_path: Path):
        """Write metadata parquet with read_id, length, etc."""
        if not HAVE_POLARS:
            # Fallback to JSON
            metadata = []
            for global_id in range(self.next_id):
                seq = self.get_sequence(global_id)
                metadata.append({
                    'read_id': global_id,
                    'length': len(seq)
                })

            with open(metadata_path.with_suffix('.json'), 'w') as f:
                json.dump(metadata, f)
            return

        # Build metadata DataFrame
        read_ids = []
        lengths = []

        for global_id in range(self.next_id):
            seq = self.get_sequence(global_id)
            read_ids.append(global_id)
            lengths.append(len(seq))

        df = pl.DataFrame({
            'read_id': read_ids,
            'length': lengths
        })

        df.write_parquet(metadata_path)

    def cleanup(self):
        """Clean up temporary files."""
        self.seq_file.close()
        if self.seq_file_path.exists():
            self.seq_file_path.unlink()


class GlobalMatrixMerger:
    """
    Merges study matrices into a global Zarr matrix.

    Process:
    1. Build global sequence vocabulary (streaming through study sequences)
    2. Create global->study ID remapping for each study
    3. Stream study matrices into global Zarr array
    """

    def __init__(self, output_dir: Path, chunk_size: int = 10000):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.chunk_size = chunk_size

        # Sequence store
        self.seq_store = MemoryMappedSequenceStore(self.output_dir / '_temp')

        # Study tracking
        self.studies: List[str] = []
        self.study_remaps: Dict[str, Dict[int, int]] = {}  # study -> {local_id: global_id}
        self.all_samples: List[str] = []
        self.study_sample_offsets: Dict[str, int] = {}  # study -> sample column offset

    def add_study(self, study_id: str, vocab_path: Path, sequences_path: Path):
        """
        Add a study's sequences to the global vocabulary.

        Returns mapping from study local IDs to global IDs.
        """
        print(f"  Adding study: {study_id}", file=sys.stderr)

        # Load study vocab
        with open(vocab_path, 'rb') as f:
            vocab_data = pickle.load(f)

        study_seq_to_id = vocab_data['seq_to_id']
        study_samples = vocab_data['samples']
        n_reads = vocab_data['n_reads']

        # Record sample offset
        self.study_sample_offsets[study_id] = len(self.all_samples)
        self.all_samples.extend(study_samples)
        self.studies.append(study_id)

        # Build local->global remapping by streaming sequences
        remap = {}

        with gzip.open(sequences_path, 'rt') as f:
            for local_id, line in enumerate(f):
                seq = line.strip()
                if not seq:
                    continue

                seq_hash = hash_sequence(seq)
                global_id, is_new = self.seq_store.add_sequence(seq, seq_hash)
                remap[local_id] = global_id

        self.study_remaps[study_id] = remap

        print(f"    {n_reads:,} study reads -> {self.seq_store.next_id:,} global reads",
              file=sys.stderr)

        return remap

    def build_global_matrix(self, study_matrices: Dict[str, Path]) -> Path:
        """
        Build the global matrix by merging all study matrices.
        """
        n_reads = self.seq_store.next_id
        n_samples = len(self.all_samples)

        print(f"\nBuilding global matrix: {n_reads:,} reads × {n_samples:,} samples",
              file=sys.stderr)

        if HAVE_ZARR:
            return self._build_zarr_matrix(study_matrices, n_reads, n_samples)
        else:
            return self._build_npz_matrix(study_matrices, n_reads, n_samples)

    def _build_zarr_matrix(self, study_matrices: Dict[str, Path],
                           n_reads: int, n_samples: int) -> Path:
        """Build chunked Zarr matrix."""
        matrix_path = self.output_dir / 'global_matrix.zarr'

        # Initialize Zarr array
        store = zarr.DirectoryStore(str(matrix_path))
        root = zarr.group(store=store, overwrite=True)

        # Create chunked array
        # Chunks: (chunk_size reads, all samples) for efficient locus queries
        counts = root.create_dataset(
            'counts',
            shape=(n_reads, n_samples),
            chunks=(self.chunk_size, n_samples),
            dtype='uint32',
            compressor=numcodecs.Blosc(cname='zstd', clevel=3, shuffle=2),
            fill_value=0
        )

        # Store metadata
        root.attrs['n_reads'] = n_reads
        root.attrs['n_samples'] = n_samples
        root.attrs['samples'] = self.all_samples
        root.attrs['chunk_size'] = self.chunk_size

        # Process each study
        for study_id, matrix_path_study in study_matrices.items():
            print(f"  Merging study: {study_id}", file=sys.stderr)

            study_matrix = sp.load_npz(matrix_path_study)
            remap = self.study_remaps[study_id]
            sample_offset = self.study_sample_offsets[study_id]

            # Convert to COO for iteration
            coo = study_matrix.tocoo()

            # Group by chunk for efficient writing
            chunk_data = {}  # chunk_idx -> list of (local_row, col, data)

            for i in range(len(coo.data)):
                local_row = coo.row[i]
                col = coo.col[i]
                val = coo.data[i]

                global_row = remap[local_row]
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
        matrix_path = self.output_dir / 'global_matrix.npz'

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

            for i in range(len(coo.data)):
                global_row = remap[coo.row[i]]
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
        with open(self.output_dir / 'samples.json', 'w') as f:
            json.dump(self.all_samples, f)

        return matrix_path

    def save_outputs(self) -> Dict[str, Path]:
        """Save all outputs and return paths."""
        self.seq_store.finalize()

        # Write FASTA
        fasta_path = self.output_dir / 'unique_reads.fasta'
        print(f"Writing FASTA: {fasta_path}", file=sys.stderr)
        self.seq_store.write_fasta(fasta_path)

        # Write metadata
        metadata_path = self.output_dir / 'read_metadata.parquet'
        print(f"Writing metadata: {metadata_path}", file=sys.stderr)
        self.seq_store.write_metadata(metadata_path)

        # Write config
        config = {
            'version': '1.0',
            'matrix_format': 'zarr' if HAVE_ZARR else 'npz',
            'chunk_size': self.chunk_size,
            'n_reads': self.seq_store.next_id,
            'n_samples': len(self.all_samples),
            'n_studies': len(self.studies),
            'studies': self.studies
        }

        config_path = self.output_dir / 'index_config.json'
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)

        # Cleanup temp files
        self.seq_store.cleanup()

        return {
            'fasta': fasta_path,
            'metadata': metadata_path,
            'config': config_path
        }


def main():
    parser = argparse.ArgumentParser(
        description='Merge study matrices into global Zarr matrix'
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

    args = parser.parse_args()

    print(f"Merging {len(args.study_dirs)} studies", file=sys.stderr)

    merger = GlobalMatrixMerger(args.output_dir, chunk_size=args.chunk_size)

    # Phase 1: Build global vocabulary
    print("\nPhase 1: Building global vocabulary", file=sys.stderr)
    study_matrices = {}

    for study_dir in args.study_dirs:
        # Find study files
        vocab_files = list(study_dir.glob('*_vocab.pkl'))
        if not vocab_files:
            print(f"Warning: No vocab file in {study_dir}, skipping", file=sys.stderr)
            continue

        study_id = vocab_files[0].stem.replace('_vocab', '')
        vocab_path = study_dir / f'{study_id}_vocab.pkl'
        seq_path = study_dir / f'{study_id}_sequences.txt.gz'
        matrix_path = study_dir / f'{study_id}_matrix.npz'

        if not all(p.exists() for p in [vocab_path, seq_path, matrix_path]):
            print(f"Warning: Incomplete study files for {study_id}, skipping", file=sys.stderr)
            continue

        merger.add_study(study_id, vocab_path, seq_path)
        study_matrices[study_id] = matrix_path

    # Phase 2: Build global matrix
    print("\nPhase 2: Building global matrix", file=sys.stderr)
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
