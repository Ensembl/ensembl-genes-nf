#!/usr/bin/env python3
"""
Convert iRibo GFF3 format to BED12.

iRibo outputs GFF3 with CDS features grouped by ID attribute.
This script groups exons by ID and creates a single BED12 entry per ORF.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
import re


def parse_gff3_attributes(attr_string):
    """Parse GFF3 attributes field into a dictionary."""
    attrs = {}
    for pair in attr_string.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs


def convert_iribo_to_bed12(input_file, output_file):
    """
    Convert iRibo GFF3 to BED12 format.

    Groups CDS features by ID and creates one BED12 entry per ORF.
    """
    # Group features by ID
    orfs = defaultdict(list)

    with open(input_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            if feature != 'CDS':
                continue

            attrs = parse_gff3_attributes(attributes)
            orf_id = attrs.get('ID', 'unknown')

            orfs[orf_id].append({
                'chrom': chrom,
                'start': int(start) - 1,  # GFF3 is 1-based, BED is 0-based
                'end': int(end),
                'strand': strand,
                'score': score if score != '.' else '0'
            })

    # Convert each ORF to BED12
    bed12_entries = []

    for orf_id, features in orfs.items():
        if not features:
            continue

        # Sort features by start position
        features = sorted(features, key=lambda x: x['start'])

        # Use first feature for shared info
        first = features[0]
        chrom = first['chrom']
        strand = first['strand']

        # Calculate transcript boundaries
        chrom_start = min(f['start'] for f in features)
        chrom_end = max(f['end'] for f in features)

        # BED12 name: use ORF ID
        name = orf_id

        # Score: use first or default
        try:
            score = int(float(first['score']))
            score = min(1000, max(0, score))
        except:
            score = 0

        # Block info
        block_count = len(features)
        block_sizes = ','.join(str(f['end'] - f['start']) for f in features)
        block_starts = ','.join(str(f['start'] - chrom_start) for f in features)

        # Use CDS bounds for thickStart/thickEnd
        thick_start = chrom_start
        thick_end = chrom_end

        # RGB color (optional, use default)
        item_rgb = '0,0,0'

        bed12_entries.append([
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
            block_sizes + ',',
            block_starts + ','
        ])

    # Sort by chr and position
    bed12_entries.sort(key=lambda x: (x[0], int(x[1])))

    # Write output
    with open(output_file, 'w') as out:
        for entry in bed12_entries:
            out.write('\t'.join(entry) + '\n')

    print(f"Converted {len(orfs)} ORFs to BED12 format", file=sys.stderr)
    print(f"Output: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Convert iRibo GFF3 to BED12 format'
    )
    parser.add_argument('input', help='Input iRibo GFF3 file')
    parser.add_argument('output', help='Output BED12 file')

    args = parser.parse_args()

    convert_iribo_to_bed12(args.input, args.output)


if __name__ == '__main__':
    main()
