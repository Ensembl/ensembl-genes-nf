#!/usr/bin/env python3
"""
Build a study-level count matrix from sample TSV files.

Takes multiple sorted TSV files (one per sample) and merges them into:
1. A sparse matrix (sequences × samples)
2. A vocabulary file (sequence hash → local ID)
3. A sequences file (for later FASTA generation)

Uses sorted two-pointer merge for memory efficiency.
Designed for 20-50 samples per study.
"""

import argparse
import gzip
import json
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterator

import numpy as np
import scipy.sparse as sp

try:
    import xxhash
    USE_XXHASH = True
except ImportError:
    import hashlib
    USE_XXHASH = False
    print("Warning: xxhash not available, using hashlib (slower)", file=sys.stderr)


def hash_sequence(seq: str) -> int:
    """Hash a sequence to a 64-bit integer."""
    if USE_XXHASH:
        return xxhash.xxh64(seq.encode()).intdigest()
    else:
        # Fallback to md5 truncated to 64 bits
        return int(hashlib.md5(seq.encode()).hexdigest()[:16], 16)


def read_sorted_tsv(tsv_path: Path) -> Iterator[Tuple[str, int]]:
    """
    Stream sorted TSV file, yielding (sequence, count) tuples.
    Assumes file is already sorted by sequence.
    """
    with open(tsv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                seq = parts[0]
                count = int(parts[1])
                yield seq, count


class StudyMatrixBuilder:
    """
    Builds a sparse count matrix for a study by merging sorted sample TSVs.

    Uses two-pointer merge to handle sorted inputs efficiently.
    Maintains sequences in sorted order for downstream global merge.
    """

    def __init__(self, study_id: str, output_dir: Path):
        self.study_id = study_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Sequence tracking
        self.sequences: List[str] = []  # Sorted list of sequences
        self.seq_to_id: Dict[int, int] = {}  # hash -> local_id

        # Count data in COO format for efficient construction
        self.rows: List[int] = []  # read_id
        self.cols: List[int] = []  # sample_idx
        self.data: List[int] = []  # count

        # Sample tracking
        self.samples: List[str] = []

    def add_sample(self, sample_id: str, tsv_path: Path):
        """
        Add a sample's reads to the study matrix.
        Uses two-pointer merge with existing sequences.
        """
        sample_idx = len(self.samples)
        self.samples.append(sample_id)

        if not self.sequences:
            # First sample - just load directly
            self._load_first_sample(sample_idx, tsv_path)
        else:
            # Merge with existing sequences
            self._merge_sample(sample_idx, tsv_path)

    def _load_first_sample(self, sample_idx: int, tsv_path: Path):
        """Load the first sample, establishing initial sequence order."""
        for seq, count in read_sorted_tsv(tsv_path):
            local_id = len(self.sequences)
            seq_hash = hash_sequence(seq)

            self.sequences.append(seq)
            self.seq_to_id[seq_hash] = local_id

            self.rows.append(local_id)
            self.cols.append(sample_idx)
            self.data.append(count)

    def _merge_sample(self, sample_idx: int, tsv_path: Path):
        """
        Merge new sample with existing sequences using two-pointer technique.
        Maintains sorted order of sequences.
        """
        # We need to insert new sequences while maintaining sort order
        # Strategy: build list of insertions, then apply at end

        new_sequences: List[Tuple[int, str]] = []  # (insert_position, sequence)
        existing_idx = 0
        n_existing = len(self.sequences)

        for new_seq, count in read_sorted_tsv(tsv_path):
            seq_hash = hash_sequence(new_seq)

            # Advance existing pointer to find position
            while existing_idx < n_existing and self.sequences[existing_idx] < new_seq:
                existing_idx += 1

            if existing_idx < n_existing and self.sequences[existing_idx] == new_seq:
                # Sequence exists - add count
                local_id = self.seq_to_id[seq_hash]
                self.rows.append(local_id)
                self.cols.append(sample_idx)
                self.data.append(count)
            else:
                # New sequence - mark for insertion
                # local_id will be assigned after all insertions
                new_sequences.append((existing_idx, new_seq, count, seq_hash))

        # Apply insertions (from end to preserve indices)
        # First, calculate final positions
        if new_sequences:
            self._apply_insertions(new_sequences, sample_idx)

    def _apply_insertions(self, new_sequences: List[Tuple], sample_idx: int):
        """
        Apply new sequence insertions, updating all data structures.
        """
        # Sort by insertion position (should already be sorted, but ensure)
        new_sequences.sort(key=lambda x: x[0])

        # Build new sequences list with insertions
        old_sequences = self.sequences
        old_seq_to_id = self.seq_to_id

        # Create mapping from old_id to new_id
        # IDs shift based on number of insertions before them
        n_insertions_before = [0] * (len(old_sequences) + 1)
        insert_positions = set(pos for pos, _, _, _ in new_sequences)

        count = 0
        for i in range(len(old_sequences) + 1):
            n_insertions_before[i] = count
            # Count insertions at this position
            count += sum(1 for pos, _, _, _ in new_sequences if pos == i)

        # Build new sequences list
        new_seq_list = []
        new_seq_to_id = {}
        insert_idx = 0

        for old_idx, seq in enumerate(old_sequences):
            # Insert any new sequences that go before this position
            while insert_idx < len(new_sequences) and new_sequences[insert_idx][0] == old_idx:
                _, ins_seq, ins_count, ins_hash = new_sequences[insert_idx]
                new_id = len(new_seq_list)
                new_seq_list.append(ins_seq)
                new_seq_to_id[ins_hash] = new_id

                # Add count entry
                self.rows.append(new_id)
                self.cols.append(sample_idx)
                self.data.append(ins_count)

                insert_idx += 1

            # Add existing sequence with new ID
            seq_hash = hash_sequence(seq)
            new_id = len(new_seq_list)
            new_seq_list.append(seq)
            new_seq_to_id[seq_hash] = new_id

        # Insert remaining new sequences at end
        while insert_idx < len(new_sequences):
            _, ins_seq, ins_count, ins_hash = new_sequences[insert_idx]
            new_id = len(new_seq_list)
            new_seq_list.append(ins_seq)
            new_seq_to_id[ins_hash] = new_id

            self.rows.append(new_id)
            self.cols.append(sample_idx)
            self.data.append(ins_count)

            insert_idx += 1

        # Remap existing row indices
        old_to_new = {}
        new_idx = 0
        insert_idx = 0
        for old_idx in range(len(old_sequences)):
            # Skip past insertions
            while insert_idx < len(new_sequences) and new_sequences[insert_idx][0] == old_idx:
                new_idx += 1
                insert_idx += 1
            old_to_new[old_idx] = new_idx
            new_idx += 1

        # Update row indices for existing data
        self.rows = [old_to_new.get(r, r) for r in self.rows]

        # Update class state
        self.sequences = new_seq_list
        self.seq_to_id = new_seq_to_id

    def save(self) -> Tuple[Path, Path, Path]:
        """
        Save study matrix and associated files.

        Returns:
            Tuple of (matrix_path, vocab_path, sequences_path)
        """
        n_reads = len(self.sequences)
        n_samples = len(self.samples)

        print(f"Study {self.study_id}: {n_reads:,} unique reads × {n_samples} samples",
              file=sys.stderr)

        # 1. Save sparse matrix (CSR format)
        matrix = sp.coo_matrix(
            (self.data, (self.rows, self.cols)),
            shape=(n_reads, n_samples),
            dtype=np.uint32
        ).tocsr()

        matrix_path = self.output_dir / f"{self.study_id}_matrix.npz"
        sp.save_npz(matrix_path, matrix)

        # 2. Save vocabulary (hash -> local_id)
        vocab_path = self.output_dir / f"{self.study_id}_vocab.pkl"
        with open(vocab_path, 'wb') as f:
            pickle.dump({
                'seq_to_id': self.seq_to_id,
                'samples': self.samples,
                'n_reads': n_reads
            }, f)

        # 3. Save sequences (gzipped, one per line, in ID order)
        sequences_path = self.output_dir / f"{self.study_id}_sequences.txt.gz"
        with gzip.open(sequences_path, 'wt') as f:
            for seq in self.sequences:
                f.write(seq + '\n')

        # 4. Save metadata
        total_counts = sum(self.data)
        metadata = {
            'study_id': self.study_id,
            'n_reads': n_reads,
            'n_samples': n_samples,
            'total_counts': total_counts,
            'samples': self.samples,
            'compression_ratio': total_counts / n_reads if n_reads > 0 else 0,
            'matrix_density': len(self.data) / (n_reads * n_samples) if n_reads * n_samples > 0 else 0
        }

        metadata_path = self.output_dir / f"{self.study_id}_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        print(f"  Matrix density: {metadata['matrix_density']:.4f}", file=sys.stderr)
        print(f"  Compression ratio: {metadata['compression_ratio']:.1f}x", file=sys.stderr)

        return matrix_path, vocab_path, sequences_path


def main():
    parser = argparse.ArgumentParser(
        description='Build study-level count matrix from sample TSVs'
    )
    parser.add_argument(
        'tsv_files',
        nargs='+',
        type=Path,
        help='Input TSV files (one per sample)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for study matrix files'
    )
    parser.add_argument(
        '-s', '--study-id',
        type=str,
        required=True,
        help='Study identifier (e.g., PRJNA12345)'
    )
    parser.add_argument(
        '--sample-ids',
        type=str,
        help='Comma-separated sample IDs (default: derived from filenames)'
    )

    args = parser.parse_args()

    # Determine sample IDs
    if args.sample_ids:
        sample_ids = args.sample_ids.split(',')
        if len(sample_ids) != len(args.tsv_files):
            print(f"Error: {len(sample_ids)} sample IDs provided but {len(args.tsv_files)} TSV files",
                  file=sys.stderr)
            sys.exit(1)
    else:
        # Derive from filenames
        sample_ids = []
        for tsv_path in args.tsv_files:
            stem = tsv_path.stem
            if stem.endswith('.tsv'):
                stem = stem[:-4]
            sample_ids.append(stem)

    print(f"Building study matrix for {args.study_id}", file=sys.stderr)
    print(f"  Samples: {len(args.tsv_files)}", file=sys.stderr)

    # Build matrix
    builder = StudyMatrixBuilder(args.study_id, args.output_dir)

    for sample_id, tsv_path in zip(sample_ids, args.tsv_files):
        print(f"  Adding sample: {sample_id}", file=sys.stderr)
        builder.add_sample(sample_id, tsv_path)

    # Save outputs
    matrix_path, vocab_path, sequences_path = builder.save()

    print(f"\nOutputs:", file=sys.stderr)
    print(f"  Matrix: {matrix_path}", file=sys.stderr)
    print(f"  Vocab: {vocab_path}", file=sys.stderr)
    print(f"  Sequences: {sequences_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
