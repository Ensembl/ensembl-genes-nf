#!/usr/bin/env python3
"""Compute gene retained fraction and append it to stable-ID SQL files."""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import Counter
from pathlib import Path


def gene_retained_fraction(decisions: list[dict[str, str]]) -> str:
    counts = Counter((row.get("type"), row.get("action")) for row in decisions)
    retained = counts[("gene", "mapped")]
    reference_total = retained + counts[("gene", "missing")]
    if reference_total == 0:
        return "n/a"
    return f"{retained / reference_total:.4f}"


def sql_identifier(name: str) -> str:
    return f"`{name}`"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decisions-tsv", type=Path, required=True)
    parser.add_argument("--sql-in", type=Path, required=True)
    parser.add_argument("--sql-out", type=Path, required=True)
    parser.add_argument("--dry-run-sql-in", type=Path, required=True)
    parser.add_argument("--dry-run-sql-out", type=Path, required=True)
    parser.add_argument("--assembly-metadata-db", required=True)
    parser.add_argument("--mapping-session-id", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.decisions_tsv.open(newline="", encoding="utf-8") as handle:
        decisions = list(csv.DictReader(handle, delimiter="\t"))

    coverage = gene_retained_fraction(decisions)

    shutil.copyfile(args.sql_in, args.sql_out)
    with args.sql_out.open("a") as handle:
        handle.write(
            f"USE {sql_identifier(args.assembly_metadata_db)};\n\n"
            "UPDATE annotation_events\n"
            "SET value = REPLACE(value, 'not_executed', "
            f"'mapped_prct:{coverage}')\n"
            f"WHERE anno_event_id = {args.mapping_session_id};\n\n"
        )

    shutil.copyfile(args.dry_run_sql_in, args.dry_run_sql_out)
    with args.dry_run_sql_out.open("a") as handle:
        handle.write(f"SELECT '{coverage}' AS gene_retained_coverage;\n")


if __name__ == "__main__":
    main()
