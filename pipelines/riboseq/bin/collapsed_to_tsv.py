#!/usr/bin/env python3
"""
Convert collapsed FASTA to TSV format, optionally partitioned by dinucleotide.

Input: Collapsed FASTA with counts in header (>read_x{count})
Output: TSV with sequence and count columns (sorted by sequence)
        When --partition is used, outputs 16 files: {sample_id}.{dinuc}.tsv

This is a simple transformation step that prepares data for study matrix building.
Partitioning at this stage means all downstream operations work within partitions,
giving ~1/16th memory usage throughout the pipeline.
"""

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List


# All possible dinucleotide prefixes (sorted for consistent ordering)
DINUCLEOTIDES = sorted([f"{a}{b}" for a in "ACGT" for b in "ACGT"])


def get_dinucleotide(seq: str) -> str:
    """Get dinucleotide prefix, defaulting to 'NN' for short/invalid sequences."""
    if len(seq) >= 2:
        prefix = seq[:2].upper()
        if prefix in DINUCLEOTIDES:
            return prefix
    return 'NN'


def parse_collapsed_fasta(fasta_path: Path) -> Dict[str, int]:
    """
    Parse collapsed FASTA file and extract sequence counts.

    Expected format:
        >read_x1234
        ATCGATCG...

    Returns dict mapping sequence -> count
    """
    reads = {}

    # Handle gzipped or plain files
    if str(fasta_path).endswith('.gz'):
        handle = gzip.open(fasta_path, 'rt')
    else:
        handle = open(fasta_path, 'r')

    with handle:
        sequence = None
        count = None

        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Parse count from header
                # Format: >read_x{count} or >anything_x{count}
                try:
                    count = int(line.split('_x')[-1])
                except (ValueError, IndexError):
                    print(f"Warning: Could not parse count from header: {line}",
                          file=sys.stderr)
                    count = 1
            else:
                # This is the sequence line
                sequence = line.upper()

                if sequence and count is not None:
                    # Accumulate counts for duplicate sequences (shouldn't happen but safe)
                    reads[sequence] = reads.get(sequence, 0) + count

                sequence = None
                count = None

    return reads


def write_sorted_tsv(reads: Dict[str, int], output_path: Path):
    """
    Write reads to TSV, sorted by sequence for efficient downstream merging.

    Output format:
        SEQUENCE\tCOUNT
    """
    # Sort by sequence for merge efficiency
    sorted_reads = sorted(reads.items(), key=lambda x: x[0])

    with open(output_path, 'w') as f:
        for sequence, count in sorted_reads:
            f.write(f"{sequence}\t{count}\n")


def write_partitioned_tsvs(reads: Dict[str, int], output_dir: Path, sample_id: str) -> List[Path]:
    """
    Write reads to 16 partition TSV files based on dinucleotide prefix.

    Output files: {sample_id}.{dinuc}.tsv (e.g., SRR123.AA.tsv, SRR123.AC.tsv, ...)

    Returns list of output paths.
    """
    # Partition reads by dinucleotide
    partitions: Dict[str, Dict[str, int]] = {di: {} for di in DINUCLEOTIDES}
    nn_partition: Dict[str, int] = {}  # For sequences that don't start with valid dinuc

    for seq, count in reads.items():
        dinuc = get_dinucleotide(seq)
        if dinuc in partitions:
            partitions[dinuc][seq] = count
        else:
            nn_partition[seq] = count

    # Warn about NN sequences (shouldn't happen with valid RNA/DNA)
    if nn_partition:
        print(f"Warning: {len(nn_partition)} sequences with invalid dinucleotide prefix (skipped)",
              file=sys.stderr)

    output_paths = []
    for dinuc in DINUCLEOTIDES:
        output_path = output_dir / f"{sample_id}.{dinuc}.tsv"
        partition_reads = partitions[dinuc]

        # Always write the file, even if empty (simplifies downstream)
        sorted_reads = sorted(partition_reads.items(), key=lambda x: x[0])
        with open(output_path, 'w') as f:
            for sequence, count in sorted_reads:
                f.write(f"{sequence}\t{count}\n")

        output_paths.append(output_path)

    return output_paths


def main():
    parser = argparse.ArgumentParser(
        description='Convert collapsed FASTA to sorted TSV(s)'
    )
    parser.add_argument(
        'input',
        type=Path,
        help='Input collapsed FASTA file (.fa or .fa.gz)'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output TSV file (default: {input_stem}.tsv). Ignored with --partition.'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('.'),
        help='Output directory for partition files (default: current directory)'
    )
    parser.add_argument(
        '--sample-id',
        type=str,
        help='Sample ID to use in output filename (default: derived from input)'
    )
    parser.add_argument(
        '--partition',
        action='store_true',
        help='Output 16 partition files by dinucleotide prefix (AA, AC, ..., TT)'
    )

    args = parser.parse_args()

    # Determine sample ID
    if args.sample_id:
        sample_id = args.sample_id
    else:
        stem = args.input.name
        if stem.endswith('.gz'):
            stem = stem[:-3]
        if stem.endswith('.fa') or stem.endswith('.fasta'):
            stem = stem.rsplit('.', 1)[0]
        sample_id = stem

    # Process
    print(f"Reading: {args.input}", file=sys.stderr)
    reads = parse_collapsed_fasta(args.input)
    print(f"Found {len(reads):,} unique sequences", file=sys.stderr)

    if args.partition:
        # Write 16 partition files
        args.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Writing 16 partition files to: {args.output_dir}", file=sys.stderr)
        output_paths = write_partitioned_tsvs(reads, args.output_dir, sample_id)

        # Report partition sizes
        partition_sizes = []
        for path in output_paths:
            with open(path, 'r') as f:
                size = sum(1 for _ in f)
            partition_sizes.append(size)

        print(f"Partition sizes: min={min(partition_sizes):,}, max={max(partition_sizes):,}, "
              f"mean={sum(partition_sizes)/len(partition_sizes):,.0f}", file=sys.stderr)
    else:
        # Write single TSV
        if args.output:
            output_path = args.output
        else:
            output_path = args.output_dir / f"{sample_id}.tsv"

        print(f"Writing: {output_path}", file=sys.stderr)
        write_sorted_tsv(reads, output_path)

    # Summary stats
    total_counts = sum(reads.values())
    print(f"Total read counts: {total_counts:,}", file=sys.stderr)
    print(f"Compression ratio: {total_counts / len(reads):.1f}x", file=sys.stderr)


if __name__ == '__main__':
    main()
