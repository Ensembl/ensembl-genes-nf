#!/usr/bin/env python3
"""
Build a study-level count matrix from sample TSV files using tournament-style merge.

Takes multiple sorted TSV files (one per sample, optionally one partition) and merges them into:
1. A sparse matrix (sequences × samples)
2. A vocabulary file (sequence hash → local ID)
3. A sequences file (for later FASTA generation)

Uses tournament-style pairwise merging for O(n log k) complexity instead of O(n²).
When used with partitioned inputs, each invocation handles one dinucleotide partition.
Designed for 20-50 samples per study.
"""

import argparse
from array import array
import gzip
import json
import pickle
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

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


VocabEntry = Union[int, List[Tuple[str, int]]]


def read_sorted_tsv(tsv_path: Path) -> Iterator[Tuple[str, int]]:
    """
    Stream sorted TSV file, yielding (sequence, count) tuples.
    Validates lexical order and coalesces duplicate adjacent sequence rows.
    """
    with open(tsv_path, 'r') as f:
        previous_seq: Optional[str] = None
        previous_count = 0

        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                raise ValueError(f"Malformed TSV line in {tsv_path}:{line_num}: {line}")

            seq = parts[0]
            try:
                count = int(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid count in {tsv_path}:{line_num}: {parts[1]}"
                ) from exc

            if count < 0:
                raise ValueError(f"Negative count in {tsv_path}:{line_num}: {count}")

            if previous_seq is not None and seq < previous_seq:
                raise ValueError(
                    f"Input TSV is not sorted by sequence at {tsv_path}:{line_num}: "
                    f"{seq!r} follows {previous_seq!r}"
                )

            if seq == previous_seq:
                previous_count += count
            else:
                if previous_seq is not None:
                    yield previous_seq, previous_count
                previous_seq = seq
                previous_count = count

        if previous_seq is not None:
            yield previous_seq, previous_count


def read_intermediate_file(path: Path) -> Iterator[Tuple[str, List[Tuple[int, int]]]]:
    """
    Read intermediate merge file.
    Format: sequence\tsample_idx:count,sample_idx:count,...

    Yields: (sequence, [(sample_idx, count), ...])
    """
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            seq = parts[0]
            counts = []
            if len(parts) > 1 and parts[1]:
                for item in parts[1].split(','):
                    sidx, cnt = item.split(':')
                    counts.append((int(sidx), int(cnt)))
            yield seq, counts


def write_intermediate_file(path: Path, data: Iterator[Tuple[str, List[Tuple[int, int]]]]):
    """
    Write intermediate merge file.
    Format: sequence\tsample_idx:count,sample_idx:count,...
    """
    with open(path, 'w') as f:
        for seq, counts in data:
            counts_str = ','.join(f"{sidx}:{cnt}" for sidx, cnt in counts)
            f.write(f"{seq}\t{counts_str}\n")


def merge_two_sources(
    source1: Iterator[Tuple[str, List[Tuple[int, int]]]],
    source2: Iterator[Tuple[str, List[Tuple[int, int]]]]
) -> Iterator[Tuple[str, List[Tuple[int, int]]]]:
    """
    Two-pointer merge of two sorted sources.
    Each source yields (sequence, [(sample_idx, count), ...])

    This is the core merge operation - O(n) for two sorted lists.
    """
    item1: Optional[Tuple[str, List[Tuple[int, int]]]] = None
    item2: Optional[Tuple[str, List[Tuple[int, int]]]] = None

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
            seq1, counts1 = item1
            seq2, counts2 = item2

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
                # Same sequence - merge counts
                merged_counts = counts1 + counts2
                yield (seq1, merged_counts)
                try:
                    item1 = next(source1)
                except StopIteration:
                    item1 = None
                try:
                    item2 = next(source2)
                except StopIteration:
                    item2 = None


def tsv_to_intermediate(tsv_path: Path, sample_idx: int) -> Iterator[Tuple[str, List[Tuple[int, int]]]]:
    """Convert a sorted TSV to intermediate format."""
    for seq, count in read_sorted_tsv(tsv_path):
        yield (seq, [(sample_idx, count)])


class TournamentMerger:
    """
    Tournament-style pairwise merger for sorted TSV files.

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

    def merge_files(self, inputs: List[Tuple[Path, int]]) -> Path:
        """
        Merge multiple (tsv_path, sample_idx) pairs using tournament merge.

        Args:
            inputs: List of (tsv_path, sample_idx) tuples

        Returns:
            Path to final merged intermediate file
        """
        if not inputs:
            raise ValueError("No inputs to merge")

        # Convert all inputs to intermediate format files
        current_files: List[Path] = []

        print(f"  Round 0: Converting {len(inputs)} TSVs to intermediate format", file=sys.stderr)

        for tsv_path, sample_idx in inputs:
            out_path = self._new_temp_file("init")
            write_intermediate_file(out_path, tsv_to_intermediate(tsv_path, sample_idx))
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


class StudyMatrixBuilder:
    """
    Builds a sparse count matrix for a study using tournament-style merge.

    All samples are merged in parallel (pairwise tournament) rather than
    sequentially, giving O(n log k) complexity.
    """

    def __init__(self, study_id: str, output_dir: Path, partition: str = None):
        self.study_id = study_id
        self.partition = partition  # e.g., "AA", "AC", etc. or None for non-partitioned
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Output file prefix
        self.prefix = f"{study_id}.{partition}" if partition else study_id

        # Will be populated after merge
        self.sequences: List[str] = []
        self.seq_to_id: Dict[int, VocabEntry] = {}
        self.samples: List[str] = []

        # COO format data. array('I') keeps large studies much smaller than
        # Python integer lists while still converting cleanly to NumPy/SciPy.
        self.rows = array('I')
        self.cols = array('I')
        self.data = array('I')

    def _add_vocab_entry(self, seq: str, local_id: int):
        """
        Record sequence -> local row mapping by hash, preserving collision buckets.

        A bare 64-bit hash is not a stable row identity. Collisions are rare, but
        if one occurs we store the colliding sequences explicitly so the global
        merge can compare the actual sequence before remapping rows.
        """
        seq_hash = hash_sequence(seq)
        existing = self.seq_to_id.get(seq_hash)

        if existing is None:
            self.seq_to_id[seq_hash] = local_id
            return

        if isinstance(existing, int):
            existing_seq = self.sequences[existing]
            if existing_seq == seq:
                return
            self.seq_to_id[seq_hash] = [(existing_seq, existing), (seq, local_id)]
            return

        for existing_seq, _ in existing:
            if existing_seq == seq:
                return
        existing.append((seq, local_id))

    @staticmethod
    def _check_uint32(value: int, label: str):
        if value > np.iinfo(np.uint32).max:
            raise OverflowError(f"{label} value exceeds uint32 range: {value}")

    def build_from_tsvs(self, sample_ids: List[str], tsv_paths: List[Path]):
        """
        Build study matrix from multiple TSV files using tournament merge.
        """
        self.samples = sample_ids
        n_samples = len(sample_ids)

        # Create temp directory for merge
        with tempfile.TemporaryDirectory(prefix=f"study_{self.prefix}_") as temp_dir:
            merger = TournamentMerger(Path(temp_dir))
            inputs = [(path, idx) for idx, path in enumerate(tsv_paths)]

            print(f"  Starting tournament merge for {n_samples} samples", file=sys.stderr)
            final_file = merger.merge_files(inputs)

            # Read final merged file and build data structures
            print(f"  Building final matrix from merged data", file=sys.stderr)
            for seq, counts in read_intermediate_file(final_file):
                local_id = len(self.sequences)

                self.sequences.append(seq)
                self._add_vocab_entry(seq, local_id)

                for sample_idx, count in counts:
                    if sample_idx < 0 or sample_idx >= n_samples:
                        raise ValueError(
                            f"Sample index out of range for {seq}: {sample_idx} "
                            f"(n_samples={n_samples})"
                        )
                    self._check_uint32(local_id, "row")
                    self._check_uint32(sample_idx, "column")
                    self._check_uint32(count, "count")
                    self.rows.append(local_id)
                    self.cols.append(sample_idx)
                    self.data.append(count)

    def save(self) -> Tuple[Path, Path, Path]:
        """
        Save study matrix and associated files.

        Returns:
            Tuple of (matrix_path, vocab_path, sequences_path)
        """
        n_reads = len(self.sequences)
        n_samples = len(self.samples)

        partition_str = f" [{self.partition}]" if self.partition else ""
        print(f"Study {self.study_id}{partition_str}: {n_reads:,} unique reads × {n_samples} samples",
              file=sys.stderr)

        # 1. Save sparse matrix (CSR format)
        row_array = np.asarray(self.rows, dtype=np.int64)
        col_array = np.asarray(self.cols, dtype=np.int64)
        data_array = np.asarray(self.data, dtype=np.uint32)

        matrix = sp.coo_matrix(
            (data_array, (row_array, col_array)),
            shape=(n_reads, n_samples),
            dtype=np.uint32
        ).tocsr()

        matrix_path = self.output_dir / f"{self.prefix}_matrix.npz"
        sp.save_npz(matrix_path, matrix)

        # 2. Save vocabulary (hash -> local_id)
        vocab_path = self.output_dir / f"{self.prefix}_vocab.pkl"
        with open(vocab_path, 'wb') as f:
            pickle.dump({
                'vocab_format_version': 2,
                'seq_to_id': self.seq_to_id,
                'samples': self.samples,
                'n_reads': n_reads,
                'partition': self.partition
            }, f)

        # 3. Save sequences (gzipped, one per line, in ID order)
        sequences_path = self.output_dir / f"{self.prefix}_sequences.txt.gz"
        with gzip.open(sequences_path, 'wt') as f:
            for seq in self.sequences:
                f.write(seq + '\n')

        # 4. Save metadata
        total_counts = sum(self.data)
        metadata = {
            'study_id': self.study_id,
            'partition': self.partition,
            'n_reads': n_reads,
            'n_samples': n_samples,
            'total_counts': total_counts,
            'samples': self.samples,
            'compression_ratio': total_counts / n_reads if n_reads > 0 else 0,
            'matrix_density': len(self.data) / (n_reads * n_samples) if n_reads * n_samples > 0 else 0
        }

        metadata_path = self.output_dir / f"{self.prefix}_metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)

        print(f"  Matrix density: {metadata['matrix_density']:.4f}", file=sys.stderr)
        print(f"  Compression ratio: {metadata['compression_ratio']:.1f}x", file=sys.stderr)

        return matrix_path, vocab_path, sequences_path


def main():
    parser = argparse.ArgumentParser(
        description='Build study-level count matrix from sample TSVs using tournament merge'
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
    parser.add_argument(
        '--partition',
        type=str,
        default=None,
        help='Partition identifier (e.g., AA, AC, ..., TT). Added to output filenames.'
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
        # Derive from filenames (remove partition suffix if present)
        sample_ids = []
        for tsv_path in args.tsv_files:
            stem = tsv_path.stem
            # Handle partitioned filenames like SRR123.AA.tsv -> SRR123
            if args.partition and stem.endswith(f".{args.partition}"):
                stem = stem[:-len(f".{args.partition}")]
            elif stem.endswith('.tsv'):
                stem = stem[:-4]
            sample_ids.append(stem)

    partition_str = f" [{args.partition}]" if args.partition else ""
    print(f"Building study matrix for {args.study_id}{partition_str}", file=sys.stderr)
    print(f"  Samples: {len(args.tsv_files)}", file=sys.stderr)

    # Build matrix using tournament merge
    builder = StudyMatrixBuilder(args.study_id, args.output_dir, partition=args.partition)
    builder.build_from_tsvs(sample_ids, args.tsv_files)

    # Save outputs
    matrix_path, vocab_path, sequences_path = builder.save()

    print(f"\nOutputs:", file=sys.stderr)
    print(f"  Matrix: {matrix_path}", file=sys.stderr)
    print(f"  Vocab: {vocab_path}", file=sys.stderr)
    print(f"  Sequences: {sequences_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
