#!/usr/bin/env python3
import argparse
from pathlib import Path
from qc_db import connect, insert_artifact, md5sum, insert_metrics

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(--db, required=True)
    ap.add_argument(--run-id, required=True)
    ap.add_argument(--sample-id, required=True)
    ap.add_argument(--report, required=True)
    ap.add_argument(--checks, required=True)
    args = ap.parse_args()

    con = connect(args.db)

    for name, path in [(getRPF_report, args.report), (getRPF_checks, args.checks)]:
        p = Path(path)
        if p.exists():
            insert_artifact(con, (
                f"{args.sample_id}:getRPF:{name}",
                args.run_id, args.sample_id, getRPF, p.name, str(p), md5sum(p), p.stat().st_size, text/plain
            ))

    # Optional simple pass fraction parse if checks is a key<TAB>value list
    try:
        ok = 0
        total = 0
        with open(args.checks) as f:
            for line in f:
                parts = line.strip().split(t)
                if len(parts) == 2:
                    v = parts[1]
                    total += 1
                    passed = 1.0 if str(v).lower() in (true, pass, 1, ok, yes) else 0.0
                    ok += passed
        if total:
            insert_metrics(con, [
                (args.run_id, args.sample_id, getRPF, getrpf.checks_pass_frac, sample, None, ok/total, None)
            ])
    except Exception:
        pass

if __name__ == __main__:
    main()
