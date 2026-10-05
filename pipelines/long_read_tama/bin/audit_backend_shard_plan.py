#!/usr/bin/env python3
"""Expand one accession's shard manifest into backend-specific planned states."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan", type=Path)
    parser.add_argument("backend")
    parser.add_argument("accession")
    parser.add_argument("skip_keys")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    skips = {value for value in args.skip_keys.split(",") if value}
    with args.plan.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"shard plan is empty: {args.plan}")
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["backend", "run_accession", "shard", "mapped_reads", "resource_class", "status"])
        for row in rows:
            shard = row.get("contig", "whole")
            key = f"{args.accession}:{shard}"
            writer.writerow([
                args.backend,
                args.accession,
                shard,
                row.get("mapped_reads", ""),
                row.get("resource_class", ""),
                "SKIPPED_BY_REQUEST" if key in skips else "EXPECTED",
            ])


if __name__ == "__main__":
    main()
