#!/usr/bin/env python3
import argparse
from pathlib import Path

from qc_db import connect, insert_artifact, insert_metrics, md5sum


STAR_KEYS = [
    ("Number of input reads", "star.total_reads"),
    ("Uniquely mapped reads number", "star.unique_mapped_reads"),
    ("Uniquely mapped reads %", "star.unique_mapped_pct"),
    ("Number of reads mapped to multiple loci", "star.multi_mapped_reads"),
    ("% of reads mapped to multiple loci", "star.multi_mapped_pct"),
]


def parse_log(path):
    vals = {}
    with open(path) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = [part.strip() for part in line.split("|", 1)]
            for star_key, metric_name in STAR_KEYS:
                if key != star_key:
                    continue
                try:
                    vals[metric_name] = float(value.replace("%", "").strip())
                except ValueError:
                    pass

    total = vals.get("star.total_reads") or 0.0
    unique = vals.get("star.unique_mapped_reads") or 0.0
    if total:
        vals["star.unique_mapped_frac"] = unique / total
    return vals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--study-id", default="unknown")
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    con = connect(args.db)
    rows = [
        (args.run_id, args.sample_id, "STAR", metric, "sample", None, float(value), None)
        for metric, value in parse_log(args.log).items()
    ]
    insert_metrics(con, rows)

    path = Path(args.log)
    insert_artifact(
        con,
        (
            f"{args.sample_id}:STAR:Log.final.out",
            args.run_id,
            args.sample_id,
            "STAR",
            path.name,
            str(path),
            md5sum(path),
            path.stat().st_size,
            "text/plain",
        ),
    )


if __name__ == "__main__":
    main()
