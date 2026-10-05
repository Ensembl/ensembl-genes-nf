#!/usr/bin/env python3
"""Resolve planned shard states against native collapse products."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("plan", type=Path)
    parser.add_argument("backend")
    parser.add_argument("accession")
    parser.add_argument("success_keys")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    successes = {value for value in args.success_keys.split(",") if value}
    with args.plan.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"shard plan is empty: {args.plan}")
    statuses = []
    for row in rows:
        shard = row.get("shard") or row.get("contig") or "whole"
        key = f"{args.backend}:{args.accession}:{shard}"
        planned = row.get("status", "EXPECTED")
        if planned == "SKIPPED_BY_REQUEST":
            status = planned
        elif key in successes:
            status = "SUCCESS"
        else:
            status = "FAILED"
        statuses.append((shard, row.get("mapped_reads", ""), row.get("resource_class", ""), status))
    successful = sum(status == "SUCCESS" for _, _, _, status in statuses)
    overall = "SUCCESS" if successful else ("EMPTY_OUTPUT" if all(status == "SKIPPED_BY_REQUEST" for *_, status in statuses) else "FAILED")
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["backend", "run_accession", "shard", "mapped_reads", "resource_class", "status", "overall_status"])
        for shard, mapped_reads, resource_class, status in statuses:
            writer.writerow([args.backend, args.accession, shard, mapped_reads, resource_class, status, overall])


if __name__ == "__main__":
    main()
