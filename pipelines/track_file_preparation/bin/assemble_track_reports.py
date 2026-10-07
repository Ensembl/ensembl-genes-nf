#!/usr/bin/env python3
"""Assemble deterministic track, entity-status, and exception reports."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows.extend(reader)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_rows(args.result)
    fields = [
        "gca_accession", "assembly_release", "entity_id", "entity_type",
        "track_type", "track_path", "source_path", "status", "sha256",
        "tool_versions", "normalization_parameters",
    ]
    rows.sort(key=lambda row: (row["gca_accession"], row["entity_id"], row["track_type"]))
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)

    statuses: dict[tuple[str, str], list[str]] = defaultdict(list)
    entities: dict[tuple[str, str], dict[str, str]] = {}
    exceptions: list[dict[str, str]] = []
    for row in rows:
        key = (row["gca_accession"], row["entity_id"])
        statuses[key].append(row["status"])
        entities[key] = row
        if row["status"] == "failed":
            exceptions.append(row)

    with Path("entity_status.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gca_accession", "assembly_release", "entity_id", "entity_type", "status"])
        for key in sorted(statuses):
            values = statuses[key]
            status = "complete" if all(value == "complete" for value in values) else (
                "partial" if any(value == "complete" for value in values) else "failed"
            )
            entity = entities[key]
            writer.writerow([key[0], entity["assembly_release"], key[1], entity["entity_type"], status])

    with Path("exceptions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in exceptions)


if __name__ == "__main__":
    main()
