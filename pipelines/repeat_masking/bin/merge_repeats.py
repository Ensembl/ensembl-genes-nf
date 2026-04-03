#!/usr/bin/env python3
"""
merge_repeats.py — Merge BED files from RepeatMasker, RED, TRF, and DUST
into a single sorted, non-redundant GFF3 and BED file.

Overlapping/adjacent intervals from different sources are merged using
bedtools merge. Source information is preserved in the GFF3 attributes.

Usage:
    merge_repeats.py --beds *.bed --out-gff3 repeats.gff3 --out-bed repeats.bed
                     [--sample SAMPLE_ID]
"""

from __future__ import annotations

import argparse
import sys
import subprocess
import tempfile
import os
from pathlib import Path
from typing import Optional


def source_from_filename(path: str) -> str:
    """Infer repeat source from filename suffix."""
    p = Path(path).name
    if '.rpt.' in p or p.endswith('.rpt.bed'):
        return 'RED'
    if '.trf.' in p:
        return 'TRF'
    if '.dust.' in p:
        return 'DUST'
    # RepeatMasker GFF — could contain type info
    return 'RepeatMasker'


def parse_bed_line(line: str, source: str) -> Optional[tuple]:
    """Parse a BED or GFF2 line, return (chrom, start, end, source, name)."""
    cols = line.rstrip('\n').split('\t')
    if len(cols) < 3:
        return None
    try:
        chrom = cols[0]
        start = int(cols[1])
        end   = int(cols[2])
        name  = cols[3] if len(cols) > 3 else source
        return chrom, start, end, source, name
    except (ValueError, IndexError):
        return None


def beds_to_sorted_bed(bed_paths: list[str], tmp_path: str) -> None:
    """Concatenate and sort all BED files into tmp_path."""
    records = []
    for path in bed_paths:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        source = source_from_filename(path)
        with open(path) as fh:
            for line in fh:
                if line.startswith('#') or line.startswith('track') or not line.strip():
                    continue
                rec = parse_bed_line(line, source)
                if rec:
                    records.append(rec)

    # Sort by chrom then start
    records.sort(key=lambda r: (r[0], r[1]))

    with open(tmp_path, 'w') as fh:
        for chrom, start, end, source, name in records:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\t{source}\n")


def merge_bed(sorted_bed: str, out_bed: str) -> None:
    """Run bedtools merge to collapse overlapping intervals."""
    cmd = [
        'bedtools', 'merge',
        '-i', sorted_bed,
        '-c', '4,6,7',    # carry name, strand, source
        '-o', 'collapse,collapse,collapse',
    ]
    with open(out_bed, 'w') as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def bed_to_gff3(bed_path: str, gff3_path: str, sample_id: str) -> int:
    """Convert merged BED to GFF3 with repeat features."""
    count = 0
    with open(bed_path) as fin, open(gff3_path, 'w') as fout:
        fout.write("##gff-version 3\n")
        for i, line in enumerate(fin, start=1):
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            chrom  = cols[0]
            start  = int(cols[1]) + 1  # GFF3 is 1-based
            end    = int(cols[2])
            name   = cols[3] if len(cols) > 3 else 'repeat'
            source = cols[6] if len(cols) > 6 else 'repeat_masking'
            feat_id = f"{sample_id}_repeat_{i:08d}"
            fout.write(
                f"{chrom}\t{source}\trepeat_region\t{start}\t{end}"
                f"\t.\t.\t.\tID={feat_id};Name={name}\n"
            )
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--beds',     nargs='+', required=True, help='Input BED files')
    parser.add_argument('--out-gff3', required=True, help='Output GFF3 file')
    parser.add_argument('--out-bed',  required=True, help='Output merged BED file')
    parser.add_argument('--sample',   default='genome', help='Sample ID for feature IDs')
    args = parser.parse_args()

    with tempfile.NamedTemporaryFile(mode='w', suffix='.sorted.bed', delete=False) as tf:
        tmp_sorted = tf.name

    try:
        print(f"[merge_repeats] Merging {len(args.beds)} BED file(s)...", flush=True)
        beds_to_sorted_bed(args.beds, tmp_sorted)

        if os.path.getsize(tmp_sorted) == 0:
            print("[merge_repeats] Warning: no repeat features found in input BEDs", flush=True)
            open(args.out_bed,  'w').close()
            open(args.out_gff3, 'w').write("##gff-version 3\n")
            return

        merge_bed(tmp_sorted, args.out_bed)
        n = bed_to_gff3(args.out_bed, args.out_gff3, args.sample)
        print(f"[merge_repeats] Wrote {n} repeat features to {args.out_gff3}", flush=True)

    finally:
        os.unlink(tmp_sorted)


if __name__ == '__main__':
    main()
