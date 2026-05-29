#!/usr/bin/env python3
import argparse
from pathlib import Path

from qc_db import connect, insert_artifact, insert_metrics, md5sum


PASS_VALUES = {"true", "pass", "passed", "1", "ok", "yes", "y"}


def parse_checks(path):
    ok = 0.0
    total = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            total += 1
            if parts[1].strip().lower() in PASS_VALUES:
                ok += 1.0
    return (ok / total) if total else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--checks", required=True)
    args = parser.parse_args()

    con = connect(args.db)

    for name, file_path in [("getRPF_report", args.report), ("getRPF_checks", args.checks)]:
        path = Path(file_path)
        if not path.exists():
            continue
        insert_artifact(
            con,
            (
                f"{args.sample_id}:getRPF:{name}",
                args.run_id,
                args.sample_id,
                "getRPF",
                path.name,
                str(path),
                md5sum(path),
                path.stat().st_size,
                "text/plain",
            ),
        )

    pass_fraction = parse_checks(args.checks)
    if pass_fraction is not None:
        insert_metrics(
            con,
            [
                (
                    args.run_id,
                    args.sample_id,
                    "getRPF",
                    "getrpf.checks_pass_frac",
                    "sample",
                    None,
                    pass_fraction,
                    None,
                )
            ],
        )


if __name__ == "__main__":
    main()
