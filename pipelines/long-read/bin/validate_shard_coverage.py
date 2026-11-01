#!/usr/bin/env python3
"""Validate that planned eligible shards conserve mapped reads from a source BAM."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--source-count", required=True, type=int)
    parser.add_argument("--shard-count", required=True, type=int)
    parser.add_argument("--shard-files", nargs="*", default=[])
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    with args.plan.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible = [row for row in rows if row.get("eligible", "").lower() == "true"]
    expected = sum(int(row["mapped_reads"]) for row in eligible)
    errors: list[str] = []
    if not eligible and args.source_count:
        errors.append("plan contains no eligible shard")
    if len(args.shard_files) != len(eligible):
        errors.append(f"expected {len(eligible)} shard files, found {len(args.shard_files)}")
    if expected != args.source_count:
        errors.append(f"plan mapped-read total {expected} differs from source {args.source_count}")
    if args.shard_count != args.source_count:
        errors.append(f"shard mapped-read total {args.shard_count} differs from source {args.source_count}")
    with args.report.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"planned_eligible_shards\t{len(eligible)}\n")
        handle.write(f"emitted_shards\t{len(args.shard_files)}\n")
        handle.write(f"source_mapped_reads\t{args.source_count}\n")
        handle.write(f"planned_mapped_reads\t{expected}\n")
        handle.write(f"shard_mapped_reads\t{args.shard_count}\n")
        handle.write(f"errors\t{len(errors)}\n")
        for error in errors:
            handle.write(f"error\t{error}\n")
    if errors:
        raise SystemExit("invalid shard coverage: " + "; ".join(errors))


if __name__ == "__main__":
    main()
