#!/usr/bin/env python3
"""
Convert collapsed FASTA to TSV format.

Input: Collapsed FASTA with counts in header (>read_x{count})
Output: TSV with sequence and count columns (sorted by sequence)

This is a simple transformation step that prepares data for study matrix building.
"""

import argparse
import gzip
import sys
from pathlib import Path


def parse_collapsed_fasta(fasta_path: Path) -> dict:
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


def write_sorted_tsv(reads: dict, output_path: Path):
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


def main():
    parser = argparse.ArgumentParser(
        description='Convert collapsed FASTA to sorted TSV'
    )
    parser.add_argument(
        'input',
        type=Path,
        help='Input collapsed FASTA file (.fa or .fa.gz)'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output TSV file (default: {input_stem}.tsv)'
    )
    parser.add_argument(
        '--sample-id',
        type=str,
        help='Sample ID to use in output filename (default: derived from input)'
    )

    args = parser.parse_args()

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        # Derive from input filename
        stem = args.input.name
        if stem.endswith('.gz'):
            stem = stem[:-3]
        if stem.endswith('.fa') or stem.endswith('.fasta'):
            stem = stem.rsplit('.', 1)[0]

        if args.sample_id:
            stem = args.sample_id

        output_path = Path(f"{stem}.tsv")

    # Process
    print(f"Reading: {args.input}", file=sys.stderr)
    reads = parse_collapsed_fasta(args.input)

    print(f"Found {len(reads):,} unique sequences", file=sys.stderr)

    print(f"Writing: {output_path}", file=sys.stderr)
    write_sorted_tsv(reads, output_path)

    # Summary stats
    total_counts = sum(reads.values())
    print(f"Total read counts: {total_counts:,}", file=sys.stderr)
    print(f"Compression ratio: {total_counts / len(reads):.1f}x", file=sys.stderr)


if __name__ == '__main__':
    main()
