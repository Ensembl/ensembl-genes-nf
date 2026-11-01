#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workload", required=True, type=Path)
    parser.add_argument("--mode", choices=["contig", "none"], required=True)
    parser.add_argument("--threshold", type=int, required=True)
    parser.add_argument("--very-large-threshold", type=int, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--plan", required=True, type=Path)
    args = parser.parse_args()
    with args.workload.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("workload is empty")
    with args.plan.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["shard_id", "contig", "reference_length", "mapped_reads", "eligible", "resource_class", "source_run_id", "shard_status"])
        if args.mode == "none":
            total = sum(int(row["mapped_reads"]) for row in rows)
            writer.writerow([f"{args.run_id}__whole_bam", "whole_bam", "", total, "true", "whole_bam", args.run_id, "PLANNED"])
            return
        for row in rows:
            mapped = int(row["mapped_reads"])
            eligible = mapped > 0
            if not eligible:
                resource = "empty"
            elif mapped >= args.very_large_threshold:
                resource = "very_large"
            elif mapped >= args.threshold:
                resource = "large"
            else:
                resource = "small"
            status = "EMPTY" if not eligible else "PLANNED"
            shard_id = f"{args.run_id}__{row['contig']}".replace("/", "_")
            writer.writerow([shard_id, row["contig"], row["reference_length"], mapped, str(eligible).lower(), resource, args.run_id, status])


if __name__ == "__main__":
    main()
