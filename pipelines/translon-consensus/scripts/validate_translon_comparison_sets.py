#!/usr/bin/env python3
"""Validate that translon DB comparison sets are present and QC-usable."""

from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path


TOOLS = ("PRICE", "RiboTIE", "ORFQuant", "iRibo", "RibORF2")
CLASSES = ("cds", "non_cds")
PANCREAS_FASTQ_SAMPLES = (
    "SRR11005875_to_79",
    "SRR11005880_to_84",
    "SRR11005885_to_89",
    "SRR11005890_to_94",
    "SRR11005895_to_99",
    "SRR11005900_to_04",
)


def read_samples(con: sqlite3.Connection) -> list[str]:
    return [
        row[0]
        for row in con.execute(
            """
            SELECT DISTINCT sample_id
            FROM translons
            WHERE run_type IN ('individual_sample', 'pooled_or_aggregate')
              AND sample_id NOT GLOB '*_fastq'
            ORDER BY sample_id
            """
        )
    ]


def load_counts(con: sqlite3.Connection) -> dict[tuple[str, str, str, str], tuple[int, int]]:
    counts: dict[tuple[str, str, str, str], tuple[int, int]] = {}
    for row in con.execute(
        """
        SELECT
            sample_id,
            source_tool,
            source_feature_class,
            input_route,
            COUNT(*) AS total_rows,
            SUM(CASE WHEN qc_status = 'pass' THEN 1 ELSE 0 END) AS pass_rows
        FROM translons
        WHERE source_feature_class IN ('cds', 'non_cds')
          AND input_route IN ('bam_to_orf', 'fastq_to_orf')
          AND sample_id NOT GLOB '*_fastq'
        GROUP BY sample_id, source_tool, source_feature_class, input_route
        """
    ):
        sample_id, tool, feature_class, route, total_rows, pass_rows = row
        counts[(sample_id, tool, feature_class, route)] = (int(total_rows or 0), int(pass_rows or 0))
    return counts


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fieldnames = list(rows[0]) if rows else ["check", "sample_id", "source_tool", "source_feature_class", "input_route", "status"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="translons.sqlite path")
    parser.add_argument("--out-dir", type=Path, required=True, help="directory for validation TSVs")
    parser.add_argument(
        "--allow-missing-fastq-tool",
        action="append",
        default=[],
        help="FASTQ-vs-BAM tool to report but not fail if missing",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    errors: list[str] = []
    rows: list[dict[str, object]] = []

    with sqlite3.connect(args.db) as con:
        samples = read_samples(con)
        counts = load_counts(con)

    allow_missing_fastq = set(args.allow_missing_fastq_tool)

    for sample_id in samples:
        for tool in TOOLS:
            for feature_class in CLASSES:
                key = (sample_id, tool, feature_class, "bam_to_orf")
                total_rows, pass_rows = counts.get(key, (0, 0))
                status = "usable_qc_pass" if pass_rows > 0 else ("present_no_qc_pass" if total_rows > 0 else "missing")
                rows.append(
                    {
                        "check": "bam_to_orf_all_samples",
                        "sample_id": sample_id,
                        "source_tool": tool,
                        "source_feature_class": feature_class,
                        "input_route": "bam_to_orf",
                        "total_rows": total_rows,
                        "pass_rows": pass_rows,
                        "status": status,
                    }
                )
                if status != "usable_qc_pass":
                    errors.append(f"BAM set not usable: {sample_id}/{tool}/{feature_class}: {status}")

    for sample_id in PANCREAS_FASTQ_SAMPLES:
        for tool in TOOLS:
            for feature_class in CLASSES:
                bam_total, bam_pass = counts.get((sample_id, tool, feature_class, "bam_to_orf"), (0, 0))
                fq_total, fq_pass = counts.get((sample_id, tool, feature_class, "fastq_to_orf"), (0, 0))
                status = "paired_qc_pass" if bam_pass > 0 and fq_pass > 0 else "not_paired_qc_pass"
                rows.append(
                    {
                        "check": "fastq_vs_bam_pancreas",
                        "sample_id": sample_id,
                        "source_tool": tool,
                        "source_feature_class": feature_class,
                        "input_route": "paired",
                        "total_rows": f"bam={bam_total};fastq={fq_total}",
                        "pass_rows": f"bam={bam_pass};fastq={fq_pass}",
                        "status": status,
                    }
                )
                if status != "paired_qc_pass" and tool not in allow_missing_fastq:
                    errors.append(f"FASTQ/BAM pair not usable: {sample_id}/{tool}/{feature_class}: bam_pass={bam_pass}, fastq_pass={fq_pass}")

    summary: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in rows:
        key = (str(row["check"]), str(row["source_tool"]), str(row["source_feature_class"]))
        rec = summary.setdefault(
            key,
            {
                "check": key[0],
                "source_tool": key[1],
                "source_feature_class": key[2],
                "usable_cells": 0,
                "total_cells": 0,
            },
        )
        rec["total_cells"] = int(rec["total_cells"]) + 1
        if row["status"] in {"usable_qc_pass", "paired_qc_pass"}:
            rec["usable_cells"] = int(rec["usable_cells"]) + 1

    detail_path = args.out_dir / "comparison_set_cell_status.tsv"
    summary_path = args.out_dir / "comparison_set_summary.tsv"
    write_tsv(detail_path, rows)
    write_tsv(summary_path, list(summary.values()))

    print(f"Wrote detail:  {detail_path}")
    print(f"Wrote summary: {summary_path}")
    if errors:
        raise SystemExit("Comparison-set validation failed:\n  - " + "\n  - ".join(errors[:50]))
    print("Comparison-set validation passed.")


if __name__ == "__main__":
    main()
