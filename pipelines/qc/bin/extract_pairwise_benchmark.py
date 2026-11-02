#!/usr/bin/env python3
"""Extract compact analysis products from a completed pairwise benchmark.

The workflow's published per-pair files are useful for audit, but are too large
to copy locally for exploratory analysis. This command streams either the
run-level TSVs or completed per-pair directories, writes small grouped
summaries, and only copies transcript-level detail for explicitly selected
pair IDs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, TextIO


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", newline=""
    )


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(path, "wt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        "w", encoding="utf-8", newline=""
    )


def rows(path: Path | list[Path]) -> Iterable[dict[str, str]]:
    paths = path if isinstance(path, list) else [path]
    for source in paths:
        with open_text(source) as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def first_row(path: Path | list[Path]) -> tuple[list[str], dict[str, str] | None]:
    source = path[0] if isinstance(path, list) else path
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], next(reader, None)


def write_tsv(path: Path, columns: list[str], records: Iterable[dict[str, str]]) -> int:
    count = 0
    with open_output(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow({column: record.get(column, "") for column in columns})
            count += 1
    return count


def integer(value: str | None) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def selected_pairs(args: argparse.Namespace) -> set[str]:
    selected = set(args.pair_id or [])
    if args.pair_list:
        selected.update(line.strip() for line in args.pair_list.read_text(encoding="utf-8").splitlines() if line.strip())
    return selected


def require(source_dir: Path, name: str) -> Path:
    path = source_dir / name
    if not path.is_file():
        raise SystemExit(f"missing required aggregate product: {path}")
    return path


def copy_selected(source: Path | list[Path], destination: Path, pair_ids: set[str]) -> int:
    columns, _ = first_row(source)
    return write_tsv(
        destination,
        columns,
        (row for row in rows(source) if not pair_ids or row.get("pair_id", "") in pair_ids),
    )


def model_inventory(source: Path) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str, str], dict[str, int | set[str]]] = {}
    for row in rows(source):
        key = (row.get("pair_id", ""), row.get("source", ""), row.get("tool", ""), row.get("technology", ""))
        group = groups.setdefault(key, {"models": 0, "loci": set(), "exons": 0, "evidence": 0})
        group["models"] += 1
        if row.get("locus_id"):
            group["loci"].add(row["locus_id"])
        group["exons"] += integer(row.get("exon_count"))
        group["evidence"] += row.get("evidence_available") == "true"
    output = []
    for (pair_id, source_name, tool, technology), group in sorted(groups.items()):
        output.append({
            "pair_id": pair_id, "source": source_name, "tool": tool, "technology": technology,
            "n_models": str(group["models"]), "n_loci": str(len(group["loci"])),
            "n_exons": str(group["exons"]), "n_models_with_evidence": str(group["evidence"]),
        })
    return output


def structural_counts(source: Path) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str, str, str], Counter[str]] = defaultdict(Counter)
    for row in rows(source):
        key = (
            row.get("pair_id", ""), row.get("tool_a", ""), row.get("tool_b", ""),
            row.get("technology_a", ""), row.get("comparison_type", ""),
        )
        groups[key][row.get("classification", "")] += 1
    output = []
    for (pair_id, tool_a, tool_b, technology_a, comparison_type), counts in sorted(groups.items()):
        for classification, count in sorted(counts.items()):
            output.append({
                "pair_id": pair_id, "tool_a": tool_a, "tool_b": tool_b,
                "technology_a": technology_a, "comparison_type": comparison_type,
                "classification": classification, "n_rows": str(count),
            })
    return output


def whole_model_counts(source: Path) -> list[dict[str, str]]:
    groups: Counter[tuple[str, str, str, str, str]] = Counter()
    for row in rows(source):
        key = (
            row.get("pair_id", ""), row.get("tool_a", ""), row.get("tool_b", ""),
            row.get("technology_a", ""), row.get("whole_model_class", ""),
        )
        groups[key] += 1
    return [
        {
            "pair_id": pair_id, "tool_a": tool_a, "tool_b": tool_b,
            "technology_a": technology_a, "whole_model_class": model_class,
            "n_models": str(count),
        }
        for (pair_id, tool_a, tool_b, technology_a, model_class), count in sorted(groups.items())
    ]


def status_counts(source: Path) -> list[dict[str, str]]:
    groups: Counter[tuple[str, str]] = Counter()
    for row in rows(source):
        groups[(row.get("metric", ""), row.get("status", ""))] += 1
    return [
        {"metric": metric, "status": status, "n_pairs": str(count)}
        for (metric, status), count in sorted(groups.items())
    ]


def aggregate_pair_summaries(source: list[Path]) -> list[dict[str, str]]:
    """Create the same comparison strata as the workflow aggregate process."""
    fields = (
        ("source_a_model_count", "total_models_a"), ("source_b_model_count", "total_models_b"),
        ("source_a_locus_count", "total_loci_a"), ("source_b_locus_count", "total_loci_b"),
        ("n_intron_comparisons", "total_intron_comparisons"),
        ("n_exact_intron_chains", "exact_intron_chains"),
        ("n_partial_intron_chains", "partial_intron_chains"),
        ("n_end_comparisons", "total_end_comparisons"),
        ("n_ends_within_threshold", "ends_within_threshold"),
        ("n_utr_comparisons", "total_utr_comparisons"),
        ("n_exact_utr_matches", "exact_utr_matches"),
        ("n_whole_model_comparisons", "total_whole_model_comparisons"),
        ("n_whole_model_matches", "whole_model_matches"),
    )
    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows(source):
        groups[(row.get("technology_a", ""), row.get("tool_a", ""), row.get("technology_b", ""), row.get("tool_b", ""))].append(row)
    output = []
    for (technology_a, tool_a, technology_b, tool_b), grouped in sorted(groups.items()):
        result = {
            "technology_a": technology_a, "tool_a": tool_a,
            "technology_b": technology_b, "tool_b": tool_b,
            "n_pairs": str(len(grouped)),
            "n_samples": str(len({row.get("sample", "") for row in grouped if row.get("sample")})),
        }
        for field, output_field in fields:
            result[output_field] = str(sum(integer(row.get(field)) for row in grouped))
        result["intron_chain_concordance"] = str(
            integer(result["exact_intron_chains"]) / integer(result["total_intron_comparisons"])
        ) if integer(result["total_intron_comparisons"]) else ""
        result["end_concordance"] = str(
            integer(result["ends_within_threshold"]) / integer(result["total_end_comparisons"])
        ) if integer(result["total_end_comparisons"]) else ""
        result["utr_concordance"] = str(
            integer(result["exact_utr_matches"]) / integer(result["total_utr_comparisons"])
        ) if integer(result["total_utr_comparisons"]) else ""
        result["whole_model_concordance"] = str(
            integer(result["whole_model_matches"]) / integer(result["total_whole_model_comparisons"])
        ) if integer(result["total_whole_model_comparisons"]) else ""
        output.append(result)
    return output


def main(args: argparse.Namespace) -> None:
    source_dir = args.results_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pair_ids = selected_pairs(args)

    aggregate_files = all((source_dir / name).is_file() for name in (
        "benchmark_manifest.tsv", "models_long.tsv", "structural_comparisons_long.tsv", "whole_model_comparisons.tsv"
    ))
    if aggregate_files:
        manifest: Path | list[Path] = require(source_dir, "benchmark_manifest.tsv")
        summary: Path | list[Path] = require(source_dir, "tool_technology_summary.tsv")
        models: Path | list[Path] = require(source_dir, "models_long.tsv")
        comparisons: Path | list[Path] = require(source_dir, "structural_comparisons_long.tsv")
        whole_models: Path | list[Path] = require(source_dir, "whole_model_comparisons.tsv")
        partial = False
    else:
        pair_dirs = sorted(path for path in source_dir.iterdir() if path.is_dir())
        manifest = sorted(path for path in source_dir.glob("*/*.summary.tsv"))
        if not pair_dirs or not manifest:
            raise SystemExit(f"no aggregate products or completed per-pair summaries found under {source_dir}")
        models = sorted(path for path in source_dir.glob("*/*.models.tsv"))
        comparisons = sorted(path for path in source_dir.glob("*/*.comparisons.tsv"))
        whole_models = sorted(path for path in source_dir.glob("*/*.whole_models.tsv"))
        if not args.summary_only and not all((models, comparisons, whole_models)):
            raise SystemExit("partial extraction needs completed model, comparison, and whole-model files; use --summary-only for pair summaries only")
        summary = []
        partial = True

    manifest_columns, _ = first_row(manifest)
    write_tsv(output_dir / "pair_summary.tsv", manifest_columns, rows(manifest))
    if partial:
        grouped_summary = aggregate_pair_summaries(manifest)
        write_tsv(
            output_dir / "tool_technology_summary.tsv",
            list(grouped_summary[0]) if grouped_summary else ["tool_a", "tool_b", "n_pairs"],
            grouped_summary,
        )
    else:
        summary_columns, _ = first_row(summary)
        write_tsv(output_dir / "tool_technology_summary.tsv", summary_columns, rows(summary))
    if not args.summary_only:
        write_tsv(
            output_dir / "model_inventory_summary.tsv",
            ["pair_id", "source", "tool", "technology", "n_models", "n_loci", "n_exons", "n_models_with_evidence"],
            model_inventory(models),
        )
        write_tsv(
            output_dir / "structural_feature_summary.tsv",
            ["pair_id", "tool_a", "tool_b", "technology_a", "comparison_type", "classification", "n_rows"],
            structural_counts(comparisons),
        )
        write_tsv(
            output_dir / "whole_model_class_counts.tsv",
            ["pair_id", "tool_a", "tool_b", "technology_a", "whole_model_class", "n_models"],
            whole_model_counts(whole_models),
        )

    status_path = source_dir / "metric_status_long.tsv"
    if status_path.is_file():
        write_tsv(
            output_dir / "metric_status_summary.tsv",
            ["metric", "status", "n_pairs"],
            status_counts(status_path),
        )

    details = output_dir / "details"
    detail_counts = {}
    if pair_ids and not args.summary_only:
        for name, source in (("models", models), ("comparisons", comparisons), ("whole_models", whole_models)):
            destination = details / f"{name}.tsv.gz"
            detail_counts[name] = copy_selected(source, destination, pair_ids)

    metadata = {
        "source_dir": str(source_dir),
        "partial_run": partial,
        "summary_only": args.summary_only,
        "selected_pair_ids": sorted(pair_ids),
        "detail_counts": detail_counts,
        "notes": [
            "Grouped files are compact summaries suitable for local plotting.",
            "Transcript-level detail is only exported when pair IDs are explicitly selected.",
            "Concordance is not truth-based recall unless an independent truth set is supplied.",
        ],
    }
    (output_dir / "bundle_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True, help="qc/pairwise_benchmark directory")
    parser.add_argument("--output-dir", type=Path, required=True, help="compact bundle directory")
    parser.add_argument("--pair-id", action="append", help="export transcript-level detail for this pair; repeatable")
    parser.add_argument("--pair-list", type=Path, help="file containing pair IDs to export")
    parser.add_argument("--summary-only", action="store_true", help="use only per-pair summary files; do not scan model/comparison tables")
    return parser


if __name__ == "__main__":
    main(parser().parse_args())
