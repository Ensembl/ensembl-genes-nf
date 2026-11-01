#!/usr/bin/env python3
"""Assemble aggregate track, status and exception reports."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = [
    "gca_accession", "assembly_release", "entity_id", "entity_type", "track_type", "track_path",
    "source_path", "status", "sha256", "tool_versions", "normalization_parameters", "aggregate_type",
    "aggregate_label", "source_entity_ids", "source_count", "aggregation_method", "metadata_json",
    "provenance_path",
]


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", nargs="*", type=Path, default=[])
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = read_rows(args.result)
    rows.sort(key=lambda row: (row.get("gca_accession", ""), row.get("entity_id", ""), row.get("track_type", "")))

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in FIELDS} for row in rows)

    expected: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    with args.plan.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            expected[(row["gca_accession"], row["aggregate_id"])].add(row["track_type"])
    observed: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        observed[(row["gca_accession"], row["entity_id"])].add(row["track_type"])

    with Path("aggregate_status.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gca_accession", "aggregate_id", "status", "expected_tracks", "observed_tracks"])
        for key in sorted(expected):
            status = "complete" if expected[key].issubset(observed[key]) else "partial"
            writer.writerow([key[0], key[1], status, ",".join(sorted(expected[key])), ",".join(sorted(observed[key]))])

    with Path("exceptions.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in FIELDS} for row in rows if row.get("status") != "complete")


if __name__ == "__main__":
    main()
