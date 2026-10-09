#!/usr/bin/env python3
"""Compare current core stable IDs with registry-derived allocations."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

from stable_id_mapping.ids import (
    exon_range_from_gene_range,
    parse_id_range,
)
from stable_id_mapping.range_check import (
    check_stable_id,
    find_duplicate_stable_ids,
    load_current_features,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--gene-range", type=parse_id_range, required=True)
    parser.add_argument("--transcript-range", type=parse_id_range, required=True)
    parser.add_argument("--translation-range", type=parse_id_range, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--invalid-tsv", type=Path, required=True)
    return parser.parse_args()


def write_invalid_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "type",
                "feature_id",
                "current_stable_id",
                "current_version",
                "status",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()

    ranges = {
        "gene": args.gene_range,
        "transcript": args.transcript_range,
        "translation": args.translation_range,
        "exon": exon_range_from_gene_range(args.gene_range),
    }

    populations = load_current_features(args.db_name)

    all_stable_ids = [
        row["stable_id"]
        for rows in populations.values()
        for row in rows
    ]
    duplicates = find_duplicate_stable_ids(all_stable_ids)

    if duplicates:
        summary = {
            "db_name": args.db_name,
            "status": "fatal_duplicates",
            "duplicate_count": len(duplicates),
            "duplicate_examples": list(duplicates[:20]),
        }
        args.output_json.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        sys.stderr.write(
            "\n*** FATAL: DUPLICATE STABLE IDs DETECTED ***\n"
            f"Database: {args.db_name}\n"
            f"Duplicate IDs: {len(duplicates)}\n"
        )
        for stable_id in duplicates[:20]:
            sys.stderr.write(f"  {stable_id}\n")

        return 2

    invalid_rows: list[dict] = []
    feature_summaries: dict[str, dict] = {}

    for feature_type, rows in populations.items():
        reason_counts: Counter[str] = Counter()

        for row in rows:
            check = check_stable_id(
                row["stable_id"],
                ranges[feature_type],
            )
            reason_counts[check.reason] += 1

            if not check.agrees:
                invalid_rows.append(
                    {
                        "type": feature_type,
                        "feature_id": row["feature_id"],
                        "current_stable_id": row["stable_id"] or "",
                        "current_version": row["version"],
                        "status": check.reason,
                    }
                )

        agreeing = reason_counts["agrees"]
        disagreeing = len(rows) - agreeing

        feature_summaries[feature_type] = {
            "total": len(rows),
            "agreeing": agreeing,
            "disagreeing": disagreeing,
            "reason_counts": dict(sorted(reason_counts.items())),
        }

    if invalid_rows:
        status = "reassignment_required"
    else:
        status = "clean"

    summary = {
        "db_name": args.db_name,
        "status": status,
        "invalid_count": len(invalid_rows),
        "features": feature_summaries,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_invalid_rows(args.invalid_tsv, invalid_rows)

    if invalid_rows:
        sys.stderr.write(
            "\n*** WARNING: STABLE-ID REASSIGNMENT IS REQUIRED ***\n"
            f"Database: {args.db_name}\n"
            f"Non-conforming IDs: {len(invalid_rows)}\n"
        )

        for feature_type, counts in feature_summaries.items():
            if counts["agreeing"] and counts["disagreeing"]:
                sys.stderr.write(
                    f"*** MIXED {feature_type.upper()} POPULATION: "
                    f"{counts['agreeing']} agreeing, "
                    f"{counts['disagreeing']} disagreeing ***\n"
                )
    else:
        print(
            f"{args.db_name}: all stable IDs agree with registry allocations"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
