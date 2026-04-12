#!/usr/bin/env python3
"""
Convert standardized ORF TSV to a minimal BED12.

Input TSV columns (tab-delimited, header required):
  sample_id,tool,chrom,start,end,strand,frame,transcript_id,orf_id,score,pval,qval,extra_json

Rules:
  - chrom: if missing/empty, fall back to transcript_id (transcript-space BED12)
  - start/end: required; rows lacking valid ints are skipped
  - name: "{orf_id}|{transcript_id}|{tool}" (use placeholders if missing)
  - score: integer 0..1000; best-effort from TSV score/pval; default 0
  - thickStart/thickEnd: full span (start..end)
  - block: single block (end-start); blockStart=0

This utility is intentionally simple; callers that can provide exon structures
or genomic mapping should emit richer BED12 upstream.
"""

import argparse
import csv
import math
from pathlib import Path


def clamp(v: int, lo: int = 0, hi: int = 1000) -> int:
    return max(lo, min(hi, v))


def score_from_fields(score_str: str | None, pval_str: str | None) -> int:
    """Map tool score/p-val to BED score (0..1000).
    Heuristics:
      - If score is numeric, linearly scale after sigmoid-like compression.
      - Else if pval is numeric in (0,1], use 1000*(1-pval).
      - Else 0.
    """
    # Try p-value first if clearly provided
    try:
        if pval_str is not None and pval_str != '' and pval_str.lower() != 'none':
            p = float(pval_str)
            if math.isfinite(p) and 0 <= p <= 1:
                return clamp(int(round(1000 * (1.0 - p))))
    except Exception:
        pass

    # Fallback to score
    try:
        if score_str is not None and score_str != '' and score_str.lower() != 'none':
            s = float(score_str)
            if math.isfinite(s):
                # squash then scale to 0..1000
                # y = 1 - exp(-|s|) in [0,1); signed not needed here
                y = 1.0 - math.exp(-abs(s))
                return clamp(int(round(1000 * y)))
    except Exception:
        pass

    return 0


def to_bed12_row(row: dict) -> str | None:
    try:
        start = int(row.get('start'))
        end = int(row.get('end'))
    except Exception:
        return None

    if end <= start:
        return None

    chrom = (row.get('chrom') or '').strip()
    tx = (row.get('transcript_id') or '').strip()
    tool = (row.get('tool') or '').strip()
    orf_id = (row.get('orf_id') or '').strip()

    if not chrom:
        # Fallback: write transcript-space BED12 using transcript_id as chrom
        chrom = tx if tx else 'unknown'

    # name field
    name_parts = [orf_id or 'orf', tx or 'tx', tool or 'tool']
    name = '|'.join(name_parts)

    score = score_from_fields(row.get('score'), row.get('pval'))
    strand = (row.get('strand') or '.').strip() or '.'

    # Full span is coding span for now
    thick_start = start
    thick_end = end

    block_count = 1
    block_sizes = f"{end - start},"
    block_starts = "0,"

    bed12 = '\t'.join([
        str(chrom),
        str(start),
        str(end),
        name,
        str(score),
        strand,
        str(thick_start),
        str(thick_end),
        '0,0,0',
        str(block_count),
        block_sizes,
        block_starts
    ])

    return bed12


def convert(tsv_path: Path, bed_path: Path) -> int:
    n_in = 0
    n_out = 0
    with open(tsv_path, newline='') as fh, open(bed_path, 'w') as out:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            n_in += 1
            bed = to_bed12_row(row)
            if bed:
                out.write(bed + '\n')
                n_out += 1
    return n_out


def main():
    ap = argparse.ArgumentParser(description='Convert standardized ORF TSV to minimal BED12')
    ap.add_argument('--input', required=True, help='Input TSV path')
    ap.add_argument('--output', required=True, help='Output BED12 path')
    args = ap.parse_args()

    tsv = Path(args.input)
    bed = Path(args.output)
    if not tsv.exists():
        raise SystemExit(f"Input TSV not found: {tsv}")
    n = convert(tsv, bed)
    # Print a tiny summary to stderr (harmless in NF logs)
    print(f"[tsv_to_bed12] Wrote {n} BED12 rows from {tsv}")


if __name__ == '__main__':
    main()

