#!/usr/bin/env python3
import argparse, json, csv
from pathlib import Path
from qc_db import connect, insert_metrics, insert_artifact, md5sum

def load_json(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return {}

def load_csv_map(p):
    out = []
    if not p or not Path(p).exists():
        return out
    with open(p) as f:
        r = csv.DictReader(f)
        for row in r:
            out.append(row)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--db', required=True)
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--json', required=False)
    ap.add_argument('--csv', required=False)
    ap.add_argument('--offsets', required=True)
    args = ap.parse_args()

    con = connect(args.db)
    rows = []

    jm = load_json(args.json) if args.json else {}
    ic = None
    for k in ('information_content','info_content'):
        if k in jm:
            ic = jm[k]
            break
    if ic is None:
        try:
            ic = jm.get('metagene',{}).get('information_content')
        except Exception:
            pass
    if ic is not None:
        try:
            rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.information_content', 'sample', None, float(ic), None))
        except Exception:
            pass
    f0 = jm.get('frame0_fraction') or (jm.get('frames',{}).get('frame0') if isinstance(jm.get('frames'), dict) else None)
    if f0 is not None:
        try:
            rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.frame0_frac_overall', 'sample', None, float(f0), None))
        except Exception:
            pass

    per_len = jm.get('per_length') or jm.get('length_metrics') or {}
    if per_len:
        for k, v in per_len.items():
            try:
                L = int(k)
            except Exception:
                continue
            reads = v.get('reads') or v.get('count') or v.get('n')
            if reads is not None:
                try:
                    rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.reads_at_length', 'length', L, float(reads), 'reads'))
                except Exception:
                    pass
            f0l = v.get('frame0_fraction') or v.get('frame0') or v.get('f0_frac')
            if f0l is not None:
                try:
                    rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.frame0_frac', 'length', L, float(f0l), None))
                except Exception:
                    pass
    else:
        for row in load_csv_map(args.csv):
            try:
                L = int(row.get('length') or row.get('read_length'))
            except Exception:
                continue
            total = row.get('total') or row.get('reads') or row.get('n')
            f0l = row.get('frame0_frac') or row.get('f0') or row.get('frame0_fraction')
            if f0l is None and 'frame0' in row and total:
                try:
                    f0l = float(row['frame0'])/float(total)
                except Exception:
                    pass
            if total is not None:
                try:
                    rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.reads_at_length', 'length', L, float(total), 'reads'))
                except Exception:
                    pass
            if f0l is not None:
                try:
                    rows.append((args.run_id, args.sample_id, 'RiboMetric', 'ribometric.frame0_frac', 'length', L, float(f0l), None))
                except Exception:
                    pass

    insert_metrics(con, rows)

    if args.json and Path(args.json).exists():
        p = Path(args.json)
        insert_artifact(con, (f'{args.sample_id}:RiboMetric:json', args.run_id, args.sample_id, 'RiboMetric', p.name, str(p), md5sum(p), p.stat().st_size, 'application/json'))
    if args.csv and Path(args.csv).exists():
        p = Path(args.csv)
        insert_artifact(con, (f'{args.sample_id}:RiboMetric:csv', args.run_id, args.sample_id, 'RiboMetric', p.name, str(p), md5sum(p), p.stat().st_size, 'text/csv'))
    p = Path(args.offsets)
    insert_artifact(con, (f'{args.sample_id}:RiboMetric:offsets', args.run_id, args.sample_id, 'RiboMetric', p.name, str(p), md5sum(p), p.stat().st_size, 'text/tab-separated-values'))

if __name__ == '__main__':
    main()
