#!/usr/bin/env python3
"""
Convert PRICE BED12 to proper BED12 format.

PRICE outputs BED12 but may have issues with block coordinates or formatting.
This script validates and fixes common issues.
"""

import argparse
import sys
from pathlib import Path


def validate_and_fix_bed12(input_file, output_file):
    """
    Validate PRICE BED12 and fix common issues.

    Fixes:
    - Trailing commas in blockSizes/blockStarts
    - Block count mismatches
    - Invalid coordinates
    """
    fixed_count = 0
    skipped_count = 0
    total_count = 0

    with open(input_file) as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('#') or line.startswith('track') or not line.strip():
                continue

            total_count += 1
            fields = line.strip().split('\t')

            if len(fields) < 12:
                print(f"Warning: Skipping line with {len(fields)} fields (expected 12)", file=sys.stderr)
                skipped_count += 1
                continue

            chrom, chrom_start, chrom_end, name, score, strand, \
                thick_start, thick_end, item_rgb, block_count, block_sizes, block_starts = fields

            try:
                # Convert to proper types
                chrom_start = int(chrom_start)
                chrom_end = int(chrom_end)
                thick_start = int(thick_start)
                thick_end = int(thick_end)
                block_count = int(block_count)

                # Clean up block_sizes and block_starts (remove trailing commas)
                block_sizes = block_sizes.rstrip(',')
                block_starts = block_starts.rstrip(',')

                # Parse blocks
                sizes = [int(x) for x in block_sizes.split(',') if x]
                starts = [int(x) for x in block_starts.split(',') if x]

                # Validate block count
                if len(sizes) != block_count or len(starts) != block_count:
                    print(f"Warning: Block count mismatch for {name}: "
                          f"declared={block_count}, sizes={len(sizes)}, starts={len(starts)}",
                          file=sys.stderr)
                    # Use actual count
                    block_count = min(len(sizes), len(starts))
                    sizes = sizes[:block_count]
                    starts = starts[:block_count]
                    fixed_count += 1

                # Validate coordinates
                if chrom_start >= chrom_end:
                    print(f"Warning: Invalid coordinates for {name}: start >= end", file=sys.stderr)
                    skipped_count += 1
                    continue

                # Validate blocks are within bounds
                valid = True
                for i in range(block_count):
                    block_abs_start = chrom_start + starts[i]
                    block_abs_end = block_abs_start + sizes[i]

                    if block_abs_end > chrom_end:
                        print(f"Warning: Block {i} extends beyond feature end for {name}", file=sys.stderr)
                        valid = False
                        break

                if not valid:
                    skipped_count += 1
                    continue

                # Ensure thickStart/thickEnd are within bounds
                thick_start = max(chrom_start, min(thick_start, chrom_end))
                thick_end = max(chrom_start, min(thick_end, chrom_end))

                # Rebuild block strings with proper formatting
                block_sizes_str = ','.join(map(str, sizes)) + ','
                block_starts_str = ','.join(map(str, starts)) + ','

                # Write fixed entry
                out.write('\t'.join([
                    chrom,
                    str(chrom_start),
                    str(chrom_end),
                    name,
                    str(score),
                    strand,
                    str(thick_start),
                    str(thick_end),
                    item_rgb,
                    str(block_count),
                    block_sizes_str,
                    block_starts_str
                ]) + '\n')

            except (ValueError, IndexError) as e:
                print(f"Warning: Error processing line: {e}", file=sys.stderr)
                skipped_count += 1
                continue

    print(f"\nProcessed {total_count} entries", file=sys.stderr)
    print(f"Fixed {fixed_count} entries", file=sys.stderr)
    print(f"Skipped {skipped_count} invalid entries", file=sys.stderr)
    print(f"Output: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Validate and fix PRICE BED12 format'
    )
    parser.add_argument('input', help='Input PRICE BED12 file')
    parser.add_argument('output', help='Output BED12 file')

    args = parser.parse_args()

    validate_and_fix_bed12(args.input, args.output)


if __name__ == '__main__':
    main()
