#!/usr/bin/env python3
"""
Filter BAM file by mapping uniqueness and junction status.

Outputs 4 BAM files based on:
- Mapping uniqueness: NH:i:1 (unique) vs NH:i:2-X (multi-mapper)
- Junction status: CIGAR contains 'N' (junction/spliced) vs no 'N' (contiguous)

Output files:
- unique_no_junction.bam: Maps to 1 genomic location, no splicing
- unique_with_junction.bam: Maps to 1 genomic location, crosses junction(s)
- multi_no_junction.bam: Maps to 2-X locations, no splicing
- multi_with_junction.bam: Maps to 2-X locations, crosses junction(s)
"""

import sys
import argparse
import pysam
from typing import Dict, Optional, Tuple
from collections import defaultdict


def has_junction(read: pysam.AlignedSegment) -> bool:
    """
    Check if read alignment crosses a junction (has 'N' in CIGAR).

    Args:
        read: pysam AlignedSegment

    Returns:
        True if CIGAR contains 'N' operation (CREF_SKIP), False otherwise
    """
    if read.cigartuples is None:
        return False
    # CIGAR operation 3 = N (reference skip / junction)
    return any(op == 3 for op, length in read.cigartuples)


def get_nh_tag(read: pysam.AlignedSegment) -> Optional[int]:
    """
    Extract NH tag (number of reported alignments) from read.

    Args:
        read: pysam AlignedSegment

    Returns:
        NH value as integer, or None if tag is missing
    """
    try:
        return read.get_tag("NH")
    except KeyError:
        return None


def get_hi_tag(read: pysam.AlignedSegment) -> Optional[int]:
    """
    Extract HI tag (hit index) from read.

    Args:
        read: pysam AlignedSegment

    Returns:
        HI value as integer, or None if tag is missing
    """
    try:
        return read.get_tag("HI")
    except KeyError:
        return None


def passes_basic_filters(read: pysam.AlignedSegment,
                         keep_primary_only: bool = True,
                         min_mapq: int = 0) -> bool:
    """
    Check if read passes basic quality filters.

    Args:
        read: pysam AlignedSegment
        keep_primary_only: If True, only keep primary alignments (HI:i:1)
        min_mapq: Minimum mapping quality threshold

    Returns:
        True if read passes filters, False otherwise
    """
    # Skip unmapped reads
    if read.is_unmapped:
        return False

    # Skip secondary alignments (0x100) and supplementary alignments (0x800)
    if read.is_secondary or read.is_supplementary:
        return False

    # Optional: keep only primary hit (HI:i:1)
    if keep_primary_only:
        hi = get_hi_tag(read)
        if hi is not None and hi != 1:
            return False

    # Optional: MAPQ filtering
    if min_mapq > 0 and read.mapping_quality < min_mapq:
        return False

    return True


def classify_read(read: pysam.AlignedSegment,
                  max_multimappers: int) -> Optional[str]:
    """
    Classify read into one of four categories based on NH tag and junction status.

    Args:
        read: pysam AlignedSegment
        max_multimappers: Maximum NH value to accept

    Returns:
        Category name as string, or None if read should be skipped
        Categories: 'unique_no_junction', 'unique_with_junction',
                   'multi_no_junction', 'multi_with_junction'
    """
    nh = get_nh_tag(read)

    # Skip reads without NH tag or exceeding threshold
    if nh is None:
        return None
    if nh > max_multimappers:
        return None

    # Check for junctions
    has_junc = has_junction(read)

    # Classify into 2x2 matrix
    if nh == 1:
        return 'unique_with_junction' if has_junc else 'unique_no_junction'
    else:  # 2 <= nh <= max_multimappers
        return 'multi_with_junction' if has_junc else 'multi_no_junction'


def open_output_bams(prefix: str, header: pysam.AlignmentHeader) -> Dict[str, pysam.AlignmentFile]:
    """
    Open 4 output BAM files for writing.

    Args:
        prefix: Output file prefix
        header: BAM header to use for output files

    Returns:
        Dictionary mapping category names to open BAM file handles
    """
    categories = [
        'unique_no_junction',
        'unique_with_junction',
        'multi_no_junction',
        'multi_with_junction'
    ]

    outputs = {}
    for category in categories:
        filename = f"{prefix}.{category}.bam"
        outputs[category] = pysam.AlignmentFile(filename, "wb", header=header)

    return outputs


def close_output_bams(outputs: Dict[str, pysam.AlignmentFile]) -> None:
    """Close all output BAM files."""
    for bam in outputs.values():
        bam.close()


def index_output_bams(prefix: str) -> None:
    """Index all output BAM files using samtools."""
    categories = [
        'unique_no_junction',
        'unique_with_junction',
        'multi_no_junction',
        'multi_with_junction'
    ]

    for category in categories:
        filename = f"{prefix}.{category}.bam"
        try:
            pysam.index(filename)
        except Exception as e:
            print(f"Warning: Failed to index {filename}: {e}", file=sys.stderr)


def filter_bam(input_bam: str,
               output_prefix: str,
               max_multimappers: int = 10,
               keep_primary_only: bool = True,
               min_mapq: int = 0,
               create_index: bool = False) -> Dict[str, int]:
    """
    Filter BAM file into 4 categories based on mapping and junction status.

    Outputs are cumulative/hierarchical:
    - unique_no_junction.bam: Most restrictive (NH=1, no junction)
    - unique_with_junction.bam: Includes above + unique with junctions
    - multi_no_junction.bam: Includes unique_no_junction + multi-mappers without junctions
    - multi_with_junction.bam: Includes all reads (least restrictive)

    Args:
        input_bam: Path to input BAM file
        output_prefix: Prefix for output files
        max_multimappers: Maximum NH value to keep (default: 10)
        keep_primary_only: Only keep primary alignments HI:i:1 (default: True)
        min_mapq: Minimum mapping quality (default: 0, disabled)
        create_index: Create BAI index files after filtering (default: False)

    Returns:
        Dictionary with read counts for each category
    """
    # Open input BAM
    inbam = pysam.AlignmentFile(input_bam, "rb")

    # Open output BAMs
    outputs = open_output_bams(output_prefix, inbam.header)

    # Statistics (count unique reads per category)
    stats = defaultdict(int)
    stats['total'] = 0

    # Stream through input BAM
    for read in inbam:
        stats['total'] += 1

        # Check basic filters
        if not passes_basic_filters(read, keep_primary_only, min_mapq):
            stats['filtered'] += 1
            continue

        # Classify read
        category = classify_read(read, max_multimappers)

        if not category:
            stats['skipped'] += 1
            continue

        # Write to appropriate output files based on cumulative hierarchy
        if category == 'unique_no_junction':
            # Most restrictive - write to all files
            outputs['unique_no_junction'].write(read)
            outputs['unique_with_junction'].write(read)
            outputs['multi_no_junction'].write(read)
            outputs['multi_with_junction'].write(read)
            stats['unique_no_junction'] += 1

        elif category == 'unique_with_junction':
            # Write to files that allow junctions
            outputs['unique_with_junction'].write(read)
            outputs['multi_with_junction'].write(read)
            stats['unique_with_junction'] += 1

        elif category == 'multi_no_junction':
            # Write to files that allow multi-mapping
            outputs['multi_no_junction'].write(read)
            outputs['multi_with_junction'].write(read)
            stats['multi_no_junction'] += 1

        elif category == 'multi_with_junction':
            # Least restrictive - only write to multi_with_junction
            outputs['multi_with_junction'].write(read)
            stats['multi_with_junction'] += 1

    # Close files
    inbam.close()
    close_output_bams(outputs)

    # Optional: create indices
    if create_index:
        print("Creating BAM indices...", file=sys.stderr)
        index_output_bams(output_prefix)

    return dict(stats)


def print_stats(stats: Dict[str, int]) -> None:
    """Print filtering statistics to stdout."""
    print("\n=== BAM Filtering Statistics ===")
    print(f"Total reads processed: {stats.get('total', 0):,}")
    print(f"\nOutput categories:")
    print(f"  Unique, no junction:   {stats.get('unique_no_junction', 0):,}")
    print(f"  Unique, with junction: {stats.get('unique_with_junction', 0):,}")
    print(f"  Multi, no junction:    {stats.get('multi_no_junction', 0):,}")
    print(f"  Multi, with junction:  {stats.get('multi_with_junction', 0):,}")
    print(f"\nFiltered out:")
    print(f"  Basic filters:         {stats.get('filtered', 0):,}")
    print(f"  Exceeds max NH:        {stats.get('skipped', 0):,}")

    total_kept = sum(stats.get(k, 0) for k in [
        'unique_no_junction', 'unique_with_junction',
        'multi_no_junction', 'multi_with_junction'
    ])
    if stats.get('total', 0) > 0:
        pct = 100 * total_kept / stats['total']
        print(f"\nTotal kept: {total_kept:,} ({pct:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        '--bam', '-b',
        required=True,
        help='Input BAM file (must be sorted)'
    )
    parser.add_argument(
        '--prefix', '-p',
        required=True,
        help='Output file prefix'
    )

    # Filtering parameters
    parser.add_argument(
        '--max-multimappers', '-m',
        type=int,
        default=10,
        help='Maximum NH value to keep (default: 10)'
    )
    parser.add_argument(
        '--min-mapq', '-q',
        type=int,
        default=0,
        help='Minimum mapping quality (default: 0, disabled)'
    )
    parser.add_argument(
        '--keep-all-hits',
        action='store_true',
        help='Keep all hits (HI:i:1,2,3...), not just primary (default: keep primary only)'
    )
    parser.add_argument(
        '--create-index',
        action='store_true',
        help='Create BAI index files after filtering'
    )

    args = parser.parse_args()

    # Run filtering
    stats = filter_bam(
        input_bam=args.bam,
        output_prefix=args.prefix,
        max_multimappers=args.max_multimappers,
        keep_primary_only=not args.keep_all_hits,
        min_mapq=args.min_mapq,
        create_index=args.create_index
    )

    # Print statistics
    print_stats(stats)


if __name__ == "__main__":
    main()
