#!/usr/bin/env python3
"""Render stable-ID SQL from decision TSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from stable_id_mapping.config import StableIdEventConfig, default_backup_prefix
from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.models import Decision
from stable_id_mapping.outputs import write_sql


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decisions-tsv", type=Path, required=True)
    parser.add_argument("--output-sql", type=Path, required=True)
    parser.add_argument("--db-name")
    parser.add_argument("--mapping-session-id", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--backup-prefix", default=default_backup_prefix())
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--replace-events-for-session", action="store_true")
    parser.add_argument("--no-translations", action="store_true")
    return parser.parse_args()


def read_decisions(path: Path) -> list[Decision]:
    decisions: list[Decision] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            decisions.append(
                Decision(
                    feature_type=row["type"],
                    action=row["action"],
                    current_stable_id=row["current_stable_id"] or None,
                    current_version=int(row["current_version"] or 0),
                    old_stable_id=row["old_stable_id"] or None,
                    old_version=int(row["old_version"] or 0),
                    new_stable_id=row["new_stable_id"] or None,
                    new_version=int(row["new_version"] or 0),
                    mapping_session_id=int(row["mapping_session_id"]),
                    score=float(row["score"] or 0),
                    reason=row["reason"],
                )
            )
    return decisions


def main() -> None:
    args = parse_args()

    decisions = read_decisions(args.decisions_tsv)
    dummy_range = parse_id_range("DUMMY:1-1")

    config = StableIdEventConfig(
        ref_gff=Path("."),
        target_gff=Path("."),
        mapped_gff=Path("."),
        report=Path("."),
        mapping_session_id=args.mapping_session_id,
        gene_range=dummy_range,
        transcript_range=dummy_range,
        translation_range=dummy_range,
        output_sql=args.output_sql,
        db_name=args.db_name,
        include_translations=not args.no_translations,
        dry_run=args.dry_run,
        backup_prefix=args.backup_prefix,
        batch_size=args.batch_size,
        replace_events_for_session=args.replace_events_for_session,
    )

    write_sql(decisions, args.output_sql, config)


if __name__ == "__main__":
    main()
