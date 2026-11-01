#!/usr/bin/env python3
"""Compare planned shard IDs with native backend audit records."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SHARD_BACKENDS = {"tama", "stringtie2", "stringtie3", "tmerge", "flair"}


def read_audits(paths: list[Path], run_id: str) -> dict[str, list[tuple[str, str]]]:
    audits: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                backend = row.get("backend", "")
                shard_id = row.get("shard_id", "")
                status = row.get("status", "")
                if row.get("run_id") and row["run_id"] != run_id:
                    continue
                if backend and shard_id:
                    audits[backend].append((shard_id, status))
    return audits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--backends", nargs="+", required=True)
    parser.add_argument("--audits", nargs="*", type=Path, default=[])
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    with args.plan.open(newline="") as handle:
        plan = list(csv.DictReader(handle, delimiter="\t"))
    if not plan:
        raise SystemExit("shard plan is empty")
    ids = [row["shard_id"] for row in plan]
    errors: list[str] = []
    if len(ids) != len(set(ids)):
        errors.append("shard IDs are not unique")
    audits = read_audits(args.audits, args.run_id)
    report_rows: list[dict[str, str]] = []
    for backend in args.backends:
        observed = defaultdict(list)
        for shard_id, status in audits.get(backend, []):
            observed[shard_id].append(status)
        if backend in SHARD_BACKENDS:
            for row in plan:
                shard_id = row["shard_id"]
                eligible = row["eligible"].lower() == "true"
                if not eligible:
                    report_rows.append({"backend": backend, "shard_id": shard_id, "planned": "SKIPPED", "attempted": "NO", "status": "SKIPPED"})
                    continue
                statuses = observed.get(shard_id, [])
                if len(statuses) == 0:
                    errors.append(f"{backend}: missing eligible shard {shard_id}")
                    attempted, status = "NO", "FAILED_MISSING"
                elif len(statuses) > 1:
                    errors.append(f"{backend}: duplicate audit records for {shard_id}")
                    attempted, status = "DUPLICATE", "FAILED_DUPLICATE"
                else:
                    attempted, status = "YES", statuses[0]
                report_rows.append({"backend": backend, "shard_id": shard_id, "planned": "PLANNED", "attempted": attempted, "status": status})
        else:
            statuses = observed.get("merged", [])
            if len(statuses) == 0:
                errors.append(f"{backend}: missing merged audit record")
                report_rows.append({"backend": backend, "shard_id": "merged", "planned": "MERGED", "attempted": "NO", "status": "FAILED_MISSING"})
            elif len(statuses) > 1:
                errors.append(f"{backend}: duplicate merged audit records")
                report_rows.append({"backend": backend, "shard_id": "merged", "planned": "MERGED", "attempted": "DUPLICATE", "status": "FAILED_DUPLICATE"})
            else:
                report_rows.append({"backend": backend, "shard_id": "merged", "planned": "MERGED", "attempted": "YES", "status": statuses[0]})
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["backend", "shard_id", "planned", "attempted", "status"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)
        handle.write(f"# errors\t{len(errors)}\n")
        for error in errors:
            handle.write(f"# error\t{error}\n")
    if errors:
        raise SystemExit("backend shard audit failed: " + "; ".join(errors))


if __name__ == "__main__":
    main()
