#!/usr/bin/env python3
"""Validate that planned eligible shards conserve mapped reads from a source BAM."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--counts", type=Path)
    parser.add_argument("--source-count", type=int)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument("--shard-files", nargs="*", default=[])
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    with args.plan.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if args.counts:
        counts = {}
        with args.counts.open() as handle:
            for line in handle:
                key, value = line.rstrip("\n").split("\t", 1)
                counts[key] = int(value)
        source_count = counts["source_mapped_reads"]
        shard_count = counts["shard_mapped_reads"]
    elif args.source_count is not None and args.shard_count is not None:
        source_count = args.source_count
        shard_count = args.shard_count
    else:
        parser.error("provide --counts or both --source-count and --shard-count")
    eligible = [row for row in rows if row.get("eligible", "").lower() == "true"]
    expected = sum(int(row["mapped_reads"]) for row in eligible)
    errors: list[str] = []
    if not eligible and source_count:
        errors.append("plan contains no eligible shard")
    if len(args.shard_files) != len(eligible):
        errors.append(f"expected {len(eligible)} shard files, found {len(args.shard_files)}")
    if expected != source_count:
        errors.append(f"plan mapped-read total {expected} differs from source {source_count}")
    if shard_count != source_count:
        errors.append(f"shard mapped-read total {shard_count} differs from source {source_count}")
    with args.report.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"planned_eligible_shards\t{len(eligible)}\n")
        handle.write(f"emitted_shards\t{len(args.shard_files)}\n")
        handle.write(f"source_mapped_reads\t{source_count}\n")
        handle.write(f"planned_mapped_reads\t{expected}\n")
        handle.write(f"shard_mapped_reads\t{shard_count}\n")
        handle.write(f"errors\t{len(errors)}\n")
        for error in errors:
            handle.write(f"error\t{error}\n")
    if errors:
        raise SystemExit("invalid shard coverage: " + "; ".join(errors))


if __name__ == "__main__":
    main()
