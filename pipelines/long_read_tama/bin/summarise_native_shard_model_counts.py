#!/usr/bin/env python3
"""Combine per-shard native model counts into one stable TSV report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = ["backend", "run_accession", "shard", "native_format", "model_count", "merge_status"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    for table in sorted(Path.cwd().rglob("shard_model_counts.tsv")):
        with table.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise SystemExit(f"unexpected columns in {table}: {reader.fieldnames}")
            rows.extend(dict(row) for row in reader)

    rows.sort(key=lambda row: (row["backend"], row["run_accession"], row["shard"]))
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
