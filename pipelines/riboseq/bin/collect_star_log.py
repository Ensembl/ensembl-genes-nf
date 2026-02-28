#!/usr/bin/env python3
import argparse
from pathlib import Path
from qc_db import connect, insert_metrics, insert_artifact, md5sum

STAR_KEYS = [
  ('Number of input reads', 'star.total_reads'),
  ('Uniquely mapped reads number', 'star.unique_mapped_reads'),
  ('Uniquely mapped reads %', 'star.unique_mapped_pct'),
  ('Number of reads mapped to multiple loci', 'star.multi_mapped_reads'),
  ('% of reads mapped to multiple loci', 'star.multi_mapped_pct'),
]

def parse_log(path):
    vals = dict()
    with open(path) as f:
        for line in f:
            if '|' not in line:
                continue
            parts = line.split('|', 1)
            k = parts[0].strip()
            v = parts[1].strip()
            for kk, outk in STAR_KEYS:
                if k == kk:
                    v_clean = v.replace('%', '').strip()
                    try:
                        vals[outk] = float(v_clean)
                    except Exception:
                        pass
    if 'star.total_reads' in vals and 'star.unique_mapped_reads' in vals:
        t = vals['star.total_reads'] or 0.0
        u = vals['star.unique_mapped_reads'] or 0.0
        vals['star.unique_mapped_frac'] = (u / t) if t else 0.0
    return vals

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', required=True)
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--study-id', default='unknown')
    ap.add_argument('--log', required=True)
    args = ap.parse_args()

    con = connect(args.db)
    metrics = parse_log(args.log)
    rows = []
    for m, v in metrics.items():
        rows.append((args.run_id, args.sample_id, 'STAR', m, 'sample', None, float(v), None))
    insert_metrics(con, rows)
    p = Path(args.log)
    insert_artifact(con, (
        '%s:STAR:Log.final.out' % (args.sample_id,),
        args.run_id, args.sample_id, 'STAR', p.name, str(p), md5sum(p), p.stat().st_size, 'text/plain'
    ))
PY}‬
