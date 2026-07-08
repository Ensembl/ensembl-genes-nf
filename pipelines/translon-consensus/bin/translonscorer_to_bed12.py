#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

def to_bed12_row(r):
    chrom = r.get('chrom') or 'chrNA'
    try:
        start = int(float(r.get('start') or 0))
        end = int(float(r.get('end') or start+1))
    except Exception:
        start, end = 0, 1
    name = r.get('orf_id') or r.get('tran_id') or 'ORF'
    try:
        score = int(float(r.get('score') or 0) * 10)
    except Exception:
        score = 0
    score = max(0, min(1000, score))
    strand = r.get('strand') or '+'
    blockCount = 1
    blockSizes = f"{max(1, end-start)},"
    blockStarts = "0,"
    thickStart = start
    thickEnd = end
    return [chrom, start, end, name, score, strand, thickStart, thickEnd, '0,0,0', blockCount, blockSizes, blockStarts]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('scored_csv', help='TranslonScorer *_orfs_scored.csv')
    ap.add_argument('out_dir', help='Output directory')
    args = ap.parse_args()

    outdir = Path(args.out_dir); outdir.mkdir(parents=True, exist_ok=True)
    src = Path(args.scored_csv)
    out = outdir / (src.stem + '.bed12')

    with open(src) as fi, open(out, 'w', newline='') as fo:
        r = csv.DictReader(fi)
        for row in r:
            fo.write('\t'.join(map(str, to_bed12_row(row)))+'\n')

    print(out)

if __name__ == '__main__':
    main()

