#!/usr/bin/env python3
"""Turn samtools idxstats output into a deterministic workload table."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--idxstats", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    with args.idxstats.open() as handle:
        for order, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] == "*":
                continue
            rows.append((order, fields[0], int(fields[1]), int(fields[2])))
    if not rows:
        raise SystemExit("idxstats contains no reference contigs")
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "reference_length", "mapped_reads", "eligible", "resource_class"])
        for _order, contig, length, mapped in rows:
            writer.writerow([contig, length, mapped, "true" if mapped else "false", "unassigned"])


if __name__ == "__main__":
    main()
