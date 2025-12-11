#!/usr/bin/env python3
"""
Convert RiboTIE GTF format to BED12.

RiboTIE outputs GTF with CDS features grouped by ORF_id.
This script groups exons by ORF_id and creates a single BED12 entry per ORF.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
import re


def parse_gtf_attributes(attr_string):
    """Parse GTF attributes field into a dictionary."""
    attrs = {}
    # Match key "value" pairs
    for match in re.finditer(r'(\w+)\s+"([^"]+)"', attr_string):
        key, value = match.groups()
        attrs[key] = value
    return attrs


def convert_ribotie_to_bed12(input_file, output_file):
    """
    Convert RiboTIE GTF to BED12 format.

    Groups CDS features by ORF_id and creates one BED12 entry per ORF.
    """
    # Group features by ORF_id
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

            attrs = parse_gtf_attributes(attributes)
            orf_id = attrs.get('ORF_id', 'unknown')

            orfs[orf_id].append({
                'chrom': chrom,
                'start': int(start) - 1,  # GTF is 1-based, BED is 0-based
                'end': int(end),
                'strand': strand,
                'gene_name': attrs.get('gene_name', ''),
                'transcript_id': attrs.get('transcript_id', ''),
                'orf_type': attrs.get('ORF_type', ''),
                'score': attrs.get('ribotie_score', '0')
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
        gene_name = first['gene_name']
        transcript_id = first['transcript_id']
        orf_type = first['orf_type']

        # Calculate transcript boundaries
        chrom_start = min(f['start'] for f in features)
        chrom_end = max(f['end'] for f in features)

        # BED12 name: combine info
        name = f"{transcript_id}_{gene_name}_{orf_type}" if gene_name else orf_id

        # Score: use average or first score
        try:
            score = int(float(first['score']) * 1000)  # Scale to 0-1000
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
            block_sizes,
            block_starts
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
        description='Convert RiboTIE GTF to BED12 format'
    )
    parser.add_argument('input', help='Input RiboTIE GTF file')
    parser.add_argument('output_dir', help='Output directory for BED12 files')

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate output filename: replace .bed with .bed12
    input_path = Path(args.input)
    if input_path.suffix == '.bed':
        output_name = input_path.stem + '.bed12'
    else:
        output_name = input_path.name + '.bed12'

    output_file = output_dir / output_name

    convert_ribotie_to_bed12(args.input, str(output_file))


if __name__ == '__main__':
    main()
