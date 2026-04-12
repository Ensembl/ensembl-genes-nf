#!/usr/bin/env python3
import argparse, json, os
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--bai', required=False, default=None)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--bam-type', required=True, choices=['transcriptome','genome'])
    ap.add_argument('--pass-lengths', required=False, default=None)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = {
        'sample_id': args.sample_id,
        'bam': str(Path(args.bam).resolve()),
        'bai': str(Path(args.bai).resolve()) if args.bai else None,
        'gtf': str(Path(args.gtf).resolve()),
        'fasta': str(Path(args.fasta).resolve()),
        'bam_type': args.bam_type,
        'pass_lengths': str(Path(args.pass_lengths).resolve()) if args.pass_lengths else None
    }
    (outdir / 'inputs.json').write_text(json.dumps(inputs, indent=2))
    (outdir / 'meta.json').write_text(json.dumps({'id': args.sample_id}, indent=2))

if __name__ == '__main__':
    main()

