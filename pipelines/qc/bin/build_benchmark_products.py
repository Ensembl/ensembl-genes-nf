#!/usr/bin/env python3
"""Create analysis-ready products from one pairwise QC comparison.

The products deliberately distinguish pairwise concordance from truth-based
recall.  A pairwise comparison can quantify agreement between two tools, but
recall requires a supplied reference or read-derived truth set.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


MODEL_COLUMNS = [
    "pair_id", "assembly", "species", "sample", "replicate", "library_protocol",
    "source", "tool", "tool_version", "technology", "source_format", "model_role",
    "transcript_id", "original_transcript_id", "original_gene_id", "locus_id",
    "gene_assignment_type", "locus_assignment_method", "model_status", "contig",
    "strand", "transcript_start", "transcript_end", "exon_count", "exon_coordinates",
    "intron_chain", "five_prime_end", "three_prime_end", "utr_coordinates",
    "cds_coordinates", "is_synthetic_locus", "read_count", "junction_support",
    "full_length_support", "evidence_source_tool", "evidence_available",
]

COMPARISON_COLUMNS = [
    "pair_id", "assembly", "species", "sample", "replicate", "library_protocol",
    "source_a", "tool_a", "tool_version_a", "technology_a", "source_b", "tool_b",
    "tool_version_b", "technology_b", "comparison_scope", "comparison_type",
    "locus_id", "feature", "source_a_transcript_id", "source_b_transcript_id",
    "classification", "shared_introns", "jaccard", "distance_bp",
    "source_a_utr_available", "source_b_utr_available", "exact_match",
]

WHOLE_MODEL_COLUMNS = [
    "pair_id", "assembly", "species", "sample", "replicate", "library_protocol",
    "source_a", "tool_a", "tool_version_a", "technology_a", "source_b", "tool_b",
    "tool_version_b", "technology_b", "comparison_scope", "locus_id",
    "source_a_transcript_id", "source_b_transcript_id", "intron_chain_exact",
    "five_prime_within_threshold", "three_prime_within_threshold", "utr_exact_match",
    "whole_model_class",
]

SUMMARY_COLUMNS = [
    "pair_id", "assembly", "species", "sample", "replicate", "library_protocol",
    "source_a", "tool_a", "tool_version_a", "technology_a", "format_a",
    "source_b", "tool_b", "tool_version_b", "technology_b", "format_b",
    "comparison_scope", "truth_source", "source_a_model_count", "source_b_model_count",
    "source_a_locus_count", "source_b_locus_count", "source_a_exon_count",
    "source_b_exon_count", "n_intron_comparisons", "n_exact_intron_chains",
    "n_partial_intron_chains", "intron_chain_concordance", "n_end_comparisons",
    "n_ends_within_threshold", "end_concordance", "n_utr_comparisons",
    "n_exact_utr_matches", "utr_concordance", "n_whole_model_comparisons",
    "n_whole_model_matches", "whole_model_concordance",
]


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def integer(value: str | None) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def rate(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.6f}" if denominator else ""


def metadata(path: str | Path) -> dict[str, str]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def common_metadata(meta: dict[str, str]) -> dict[str, str]:
    return {
        "pair_id": meta.get("pair_id", ""),
        "assembly": meta.get("assembly", ""),
        "species": meta.get("species", ""),
        "sample": meta.get("sample", ""),
        "replicate": meta.get("replicate", ""),
        "library_protocol": meta.get("library_protocol", ""),
        "source_a": meta.get("source_a", ""),
        "tool_a": meta.get("tool_a", meta.get("source_a", "")),
        "tool_version_a": meta.get("tool_version_a", ""),
        "technology_a": meta.get("technology_a", ""),
        "source_b": meta.get("source_b", ""),
        "tool_b": meta.get("tool_b", meta.get("source_b", "")),
        "tool_version_b": meta.get("tool_version_b", ""),
        "technology_b": meta.get("technology_b", ""),
    }


def model_products(meta: dict[str, str], source_a_path: str, source_b_path: str,
                   evidence_a_path: str = "", evidence_b_path: str = "") -> list[dict[str, str]]:
    base = common_metadata(meta)
    rows = []
    evidence_by_source = {}
    for source, path in ((base["source_a"], evidence_a_path), (base["source_b"], evidence_b_path)):
        evidence_by_source[source] = {}
        if path and Path(path).is_file() and Path(path).stat().st_size:
            evidence_by_source[source] = {row.get("transcript_id", ""): row for row in read_tsv(path)}
    for source, role, tool, version, technology, path in (
        (base["source_a"], "source_a", base["tool_a"], base["tool_version_a"], base["technology_a"], source_a_path),
        (base["source_b"], "source_b", base["tool_b"], base["tool_version_b"], base["technology_b"], source_b_path),
    ):
        for row in read_tsv(path):
            evidence = evidence_by_source.get(source, {}).get(row.get("transcript_id", ""), {})
            output = {
                **base,
                "source": source,
                "tool": tool,
                "tool_version": version,
                "technology": technology,
                "source_format": row.get("source_format", ""),
                "model_role": role,
                "is_synthetic_locus": str(row.get("gene_assignment_type", "") == "synthetic_locus").lower(),
                "read_count": evidence.get("read_count", ""),
                "junction_support": evidence.get("junction_support", ""),
                "full_length_support": evidence.get("full_length_support", ""),
                "evidence_source_tool": evidence.get("source_tool", ""),
                "evidence_available": str(bool(evidence)).lower(),
            }
            output.update({column: row.get(column, "") for column in MODEL_COLUMNS if column in row})
            rows.append(output)
    return rows


def comparison_products(meta: dict[str, str], intron_path: str, ends_path: str, utr_path: str) -> list[dict[str, str]]:
    base = common_metadata(meta)
    rows: list[dict[str, str]] = []

    for row in read_tsv(intron_path):
        rows.append({
            **base, "comparison_scope": "pairwise_concordance", "comparison_type": "intron_chain",
            "locus_id": row.get("locus_id", ""), "feature": "intron_chain",
            "source_a_transcript_id": row.get("source_a_transcript_id", ""),
            "source_b_transcript_id": row.get("source_b_transcript_id", ""),
            "classification": row.get("classification", ""),
            "shared_introns": row.get("shared_introns", ""), "jaccard": row.get("jaccard", ""),
        })
    for row in read_tsv(ends_path):
        rows.append({
            **base, "comparison_scope": "pairwise_concordance", "comparison_type": "transcript_end",
            "locus_id": row.get("locus_id", ""), "feature": row.get("end", ""),
            "source_a_transcript_id": row.get("source_a_transcript_id", ""),
            "source_b_transcript_id": row.get("source_b_transcript_id", ""),
            "classification": row.get("classification", ""), "distance_bp": row.get("distance_bp", ""),
        })
    for row in read_tsv(utr_path):
        rows.append({
            **base, "comparison_scope": "pairwise_concordance", "comparison_type": "utr_structure",
            "locus_id": row.get("locus_id", ""), "feature": "utr",
            "source_a_transcript_id": row.get("source_a_transcript_id", ""),
            "source_b_transcript_id": row.get("source_b_transcript_id", ""),
            "classification": row.get("classification", ""),
            "source_a_utr_available": row.get("source_a_utr_available", ""),
            "source_b_utr_available": row.get("source_b_utr_available", ""),
            "exact_match": row.get("exact_match", ""),
        })
    return rows


def whole_model_products(meta: dict[str, str], intron_path: str, ends_path: str, utr_path: str) -> list[dict[str, str]]:
    base = common_metadata(meta)
    end_matches: defaultdict[tuple[str, str, str], dict[str, bool]] = defaultdict(dict)
    for row in read_tsv(ends_path):
        key = (row.get("locus_id", ""), row.get("source_a_transcript_id", ""), row.get("source_b_transcript_id", ""))
        end_matches[key][row.get("end", "")] = row.get("classification") == "within_threshold"
    utr_matches = {}
    for row in read_tsv(utr_path):
        key = (row.get("locus_id", ""), row.get("source_a_transcript_id", ""), row.get("source_b_transcript_id", ""))
        utr_matches[key] = row

    rows = []
    for row in read_tsv(intron_path):
        key = (row.get("locus_id", ""), row.get("source_a_transcript_id", ""), row.get("source_b_transcript_id", ""))
        end_status = end_matches.get(key, {})
        five_prime = end_status.get("five_prime_end", False)
        three_prime = end_status.get("three_prime_end", False)
        utr = utr_matches.get(key, {})
        utr_exact = utr.get("exact_match", "")
        intron_exact = row.get("classification") == "exact"
        has_partner = bool(row.get("source_b_transcript_id"))
        if not has_partner:
            model_class = "missing_partner"
        elif intron_exact and five_prime and three_prime and utr_exact == "1":
            model_class = "whole_model_with_utr"
        elif intron_exact and five_prime and three_prime:
            model_class = "whole_model_structure_only"
        elif intron_exact:
            model_class = "intron_chain_only"
        else:
            model_class = "partial_or_discordant"
        rows.append({
            **base, "comparison_scope": "pairwise_concordance", "locus_id": row.get("locus_id", ""),
            "source_a_transcript_id": row.get("source_a_transcript_id", ""),
            "source_b_transcript_id": row.get("source_b_transcript_id", ""),
            "intron_chain_exact": str(int(intron_exact)),
            "five_prime_within_threshold": str(int(five_prime)),
            "three_prime_within_threshold": str(int(three_prime)),
            "utr_exact_match": utr_exact,
            "whole_model_class": model_class,
        })
    return rows


def summary_product(meta: dict[str, str], models: list[dict[str, str]], structural_path: str,
                    whole_models: list[dict[str, str]]) -> dict[str, str]:
    base = common_metadata(meta)
    by_source: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in models:
        by_source[row.get("source", "")].append(row)

    def source_count(source: str, field: str | None = None) -> int:
        rows = by_source.get(source, [])
        if field is None:
            return len(rows)
        return sum(integer(row.get(field)) for row in rows)

    def locus_count(source: str) -> int:
        return len({row.get("locus_id") for row in by_source.get(source, []) if row.get("locus_id")})

    structural = read_tsv(structural_path)
    structural = structural[0] if structural else {}
    n_intron = integer(structural.get("n_intron_comparisons"))
    exact_intron = integer(structural.get("n_exact_intron_chains"))
    partial_intron = integer(structural.get("n_partial_intron_chains"))
    n_ends = integer(structural.get("n_end_comparisons"))
    within_ends = integer(structural.get("n_ends_within_threshold"))
    n_utr = integer(structural.get("n_utr_comparisons"))
    exact_utr = integer(structural.get("n_exact_utr_matches"))
    whole_model_matches = sum(row.get("whole_model_class") in {"whole_model_with_utr", "whole_model_structure_only"} for row in whole_models)

    return {
        **base,
        "format_a": meta.get("format_a", ""), "format_b": meta.get("format_b", ""),
        "comparison_scope": "pairwise_concordance", "truth_source": "none",
        "source_a_model_count": str(source_count(base["source_a"])),
        "source_b_model_count": str(source_count(base["source_b"])),
        "source_a_locus_count": str(locus_count(base["source_a"])),
        "source_b_locus_count": str(locus_count(base["source_b"])),
        "source_a_exon_count": str(source_count(base["source_a"], "exon_count")),
        "source_b_exon_count": str(source_count(base["source_b"], "exon_count")),
        "n_intron_comparisons": str(n_intron), "n_exact_intron_chains": str(exact_intron),
        "n_partial_intron_chains": str(partial_intron),
        "intron_chain_concordance": rate(exact_intron, n_intron),
        "n_end_comparisons": str(n_ends), "n_ends_within_threshold": str(within_ends),
        "end_concordance": rate(within_ends, n_ends),
        "n_utr_comparisons": str(n_utr), "n_exact_utr_matches": str(exact_utr),
        "utr_concordance": rate(exact_utr, n_utr),
        "n_whole_model_comparisons": str(len(whole_models)),
        "n_whole_model_matches": str(whole_model_matches),
        "whole_model_concordance": rate(whole_model_matches, len(whole_models)),
    }


def build_pair(args: argparse.Namespace) -> None:
    meta = metadata(args.metadata)
    models = model_products(meta, args.source_a, args.source_b, args.evidence_a, args.evidence_b)
    comparisons = comparison_products(meta, args.intron, args.ends, args.utr)
    whole_models = whole_model_products(meta, args.intron, args.ends, args.utr)
    summary = summary_product(meta, models, args.structural, whole_models)
    write_tsv(args.output_models, MODEL_COLUMNS, models)
    write_tsv(args.output_comparisons, COMPARISON_COLUMNS, comparisons)
    write_tsv(args.output_whole_models, WHOLE_MODEL_COLUMNS, whole_models)
    write_tsv(args.output_summary, SUMMARY_COLUMNS, [summary])


def merge_tsv(paths: list[str], output: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    columns: list[str] = []
    for path in paths:
        file_rows = read_tsv(path)
        if not columns and file_rows:
            columns = list(file_rows[0])
        rows.extend(file_rows)
    write_tsv(output, columns, rows)
    return rows


def markdown(value: str | int | None) -> str:
    return str(value or "").replace("|", "\\|").replace("\n", " ")


def benchmark_report(
    models: list[dict[str, str]],
    whole_models: list[dict[str, str]],
    summaries: list[dict[str, str]],
    grouped: list[dict[str, str]],
    statuses: list[dict[str, str]],
) -> str:
    """Render a conservative run-level report from the aggregate products."""
    tools = sorted({row.get("tool", "") for row in models if row.get("tool")})
    technologies = sorted({row.get("technology", "") for row in models if row.get("technology")})
    samples = sorted({row.get("sample", "") for row in summaries if row.get("sample")})
    assignment_counts: defaultdict[str, int] = defaultdict(int)
    for row in models:
        assignment_counts[row.get("gene_assignment_type", "unassigned")] += 1
    class_counts: defaultdict[str, int] = defaultdict(int)
    for row in whole_models:
        class_counts[row.get("whole_model_class", "")] += 1
    evidence_models = sum(row.get("evidence_available") == "true" for row in models)
    status_counts: defaultdict[str, int] = defaultdict(int)
    for row in statuses:
        status_counts[row.get("status", "")] += 1

    lines = [
        "# Long-read transcript reconstruction benchmark report",
        "",
        "## Executive summary",
        "",
        f"This run contains **{len(summaries)} completed pairwise comparisons**, "
        f"{len(models)} canonical transcript models, and {len(whole_models)} "
        "whole-model comparison records.",
        "",
        "The primary results are **pairwise concordance** between reconstruction "
        "tools. They show how model inventories and structures differ, but they "
        "do not establish recall, precision, F1, biological correctness, or a "
        "best tool without a declared truth set.",
        "",
        "## Scope and provenance",
        "",
        f"- Tools represented: {', '.join(markdown(value) for value in tools) or 'not declared'}",
        f"- Technologies represented: {', '.join(markdown(value) for value in technologies) or 'not declared'}",
        f"- Samples represented: {', '.join(markdown(value) for value in samples) or 'not declared'}",
        "- Loci: constructed from the union of both models in each pair, preserving the same assignment for both sources.",
        "- Model identity: original transcript IDs are retained in the canonical and benchmark tables.",
        "",
        "## Tool and technology comparisons",
        "",
        "The following rows are the directly analysable comparison strata. Counts are summed across completed pairs; concordance is weighted by the number of comparisons.",
        "",
        "| Technology A | Tool A | Technology B | Tool B | Pairs | Models A | Models B | Intron concordance | End concordance | UTR concordance | Whole-model concordance |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in grouped:
        lines.append(
            "| " + " | ".join([
                markdown(row.get("technology_a")), markdown(row.get("tool_a")),
                markdown(row.get("technology_b")), markdown(row.get("tool_b")),
                markdown(row.get("n_pairs")), markdown(row.get("total_models_a")),
                markdown(row.get("total_models_b")),
                markdown(row.get("intron_chain_concordance")),
                markdown(row.get("end_concordance")),
                markdown(row.get("utr_concordance")),
                markdown(row.get("whole_model_concordance")),
            ]) + " |"
        )
    if not grouped:
        lines.append("| No completed comparison strata | | | | 0 | 0 | 0 | | | | |")

    lines.extend([
        "",
        "## Model production and locus assignment",
        "",
        "The model tables support model-count, exon-count, locus-count, transcript-structure, and evidence-support analyses. Synthetic loci are labels for structural grouping and are not biological gene annotations.",
        "",
        "| Assignment type | Canonical models |",
        "|---|---:|",
    ])
    for assignment, count in sorted(assignment_counts.items()):
        lines.append(f"| {markdown(assignment)} | {count} |")
    lines.extend([
        "",
        f"Models with a transcript-level evidence join: **{evidence_models} / {len(models)}**.",
        "Evidence fields are only interpreted when joined by canonical transcript ID; run-level totals are not treated as transcript evidence.",
        "",
        "## Whole-model classes",
        "",
        "Whole-model classes combine intron-chain and transcript-end agreement, with UTR-complete agreement reported separately where UTR fields exist.",
        "",
        "| Class | Models |",
        "|---|---:|",
    ])
    for model_class, count in sorted(class_counts.items()):
        lines.append(f"| {markdown(model_class)} | {count} |")
    if not class_counts:
        lines.append("| No completed whole-model comparisons | 0 |")

    lines.extend([
        "",
        "## Metric capability and run status",
        "",
        "The specification requires explicit distinction between enabled, not applicable, missing input, and failed metrics. The aggregate status counts below summarise the per-pair `metric_status.tsv` files.",
        "",
        "| Status | Entries |",
        "|---|---:|",
    ])
    for status, count in sorted(status_counts.items()):
        lines.append(f"| {markdown(status)} | {count} |")
    if not status_counts:
        lines.append("| No metric-status files supplied | 0 |")

    lines.extend([
        "",
        "The report follows the capability policy: transcript/exon/locus counts, intron-chain comparison, transcript-end comparison, and structural UTR comparison are valid for canonical and synthetic-locus inputs. Biological gene presence, gene pairing, CDS/coding integrity, protein support, and NMD claims require the corresponding biological inputs and are not inferred from BED12 or synthetic annotations.",
        "",
        "## Interpretation limits and required next analysis",
        "",
        "1. Pairwise concordance is not recall. To estimate intron-chain, whole-model, UTR, or boundary recall, supply a declared read-derived truth table or independent reference truth set.",
        "2. Tool comparisons must be stratified by technology, sample, depth/coverage, and protocol before attributing differences to the reconstruction algorithm.",
        "3. Model count differences should be reported alongside structural agreement and evidence support; model count alone is not a quality measure.",
        "4. Any missing or failed metric-status entries must be resolved or disclosed before using the results for tool-selection decisions.",
        "",
        "## Reproducible products",
        "",
        "- `models_long.tsv`: canonical model inventory with provenance and optional transcript-level evidence.",
        "- `structural_comparisons_long.tsv`: intron-chain, transcript-end, and UTR comparisons.",
        "- `whole_model_comparisons.tsv`: combined structural model classes.",
        "- `benchmark_manifest.tsv`: one summary row per completed pair.",
        "- `tool_technology_summary.tsv`: comparison strata used in this report.",
        "- `metric_status_long.tsv`: per-pair metric applicability and failures.",
        "- `metric_definitions.tsv`: metric denominators and truth-set requirements.",
        "",
    ])
    return "\n".join(lines)


def aggregate(args: argparse.Namespace) -> None:
    models = merge_tsv(args.models, "models_long.tsv")
    comparisons = merge_tsv(args.comparisons, "structural_comparisons_long.tsv")
    whole_models = merge_tsv(args.whole_models, "whole_model_comparisons.tsv")
    summaries = merge_tsv(args.summaries, "benchmark_manifest.tsv")
    statuses = merge_tsv(args.status, "metric_status_long.tsv")
    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in summaries:
        groups[(row.get("technology_a", ""), row.get("tool_a", ""), row.get("technology_b", ""), row.get("tool_b", ""))].append(row)

    summary_rows = []
    for (technology_a, tool_a, technology_b, tool_b), rows in sorted(groups.items()):
        summary_rows.append({
            "technology_a": technology_a, "tool_a": tool_a,
            "technology_b": technology_b, "tool_b": tool_b,
            "n_pairs": str(len(rows)),
            "n_samples": str(len({row.get("sample", "") for row in rows if row.get("sample")})),
            "total_models_a": str(sum(integer(row.get("source_a_model_count")) for row in rows)),
            "total_models_b": str(sum(integer(row.get("source_b_model_count")) for row in rows)),
            "total_loci_a": str(sum(integer(row.get("source_a_locus_count")) for row in rows)),
            "total_loci_b": str(sum(integer(row.get("source_b_locus_count")) for row in rows)),
            "total_intron_comparisons": str(sum(integer(row.get("n_intron_comparisons")) for row in rows)),
            "exact_intron_chains": str(sum(integer(row.get("n_exact_intron_chains")) for row in rows)),
            "partial_intron_chains": str(sum(integer(row.get("n_partial_intron_chains")) for row in rows)),
            "intron_chain_concordance": rate(sum(integer(row.get("n_exact_intron_chains")) for row in rows), sum(integer(row.get("n_intron_comparisons")) for row in rows)),
            "total_end_comparisons": str(sum(integer(row.get("n_end_comparisons")) for row in rows)),
            "ends_within_threshold": str(sum(integer(row.get("n_ends_within_threshold")) for row in rows)),
            "end_concordance": rate(sum(integer(row.get("n_ends_within_threshold")) for row in rows), sum(integer(row.get("n_end_comparisons")) for row in rows)),
            "total_utr_comparisons": str(sum(integer(row.get("n_utr_comparisons")) for row in rows)),
            "exact_utr_matches": str(sum(integer(row.get("n_exact_utr_matches")) for row in rows)),
            "utr_concordance": rate(sum(integer(row.get("n_exact_utr_matches")) for row in rows), sum(integer(row.get("n_utr_comparisons")) for row in rows)),
            "total_whole_model_comparisons": str(sum(integer(row.get("n_whole_model_comparisons")) for row in rows)),
            "whole_model_matches": str(sum(integer(row.get("n_whole_model_matches")) for row in rows)),
            "whole_model_concordance": rate(sum(integer(row.get("n_whole_model_matches")) for row in rows), sum(integer(row.get("n_whole_model_comparisons")) for row in rows)),
        })
    summary_columns = list(summary_rows[0]) if summary_rows else ["technology_a", "tool_a", "technology_b", "tool_b", "n_pairs"]
    write_tsv("tool_technology_summary.tsv", summary_columns, summary_rows)
    write_tsv("metric_definitions.tsv", ["metric", "scope", "numerator", "denominator", "truth_required", "description"], [
        {"metric": "intron_chain_concordance", "scope": "pairwise_concordance", "numerator": "exact intron-chain matches", "denominator": "intron comparisons", "truth_required": "yes for recall", "description": "Agreement of source A models with their best source B partner."},
        {"metric": "end_concordance", "scope": "pairwise_concordance", "numerator": "ends within configured threshold", "denominator": "end comparisons", "truth_required": "yes for recall", "description": "Agreement of transcript boundaries between paired models."},
        {"metric": "utr_concordance", "scope": "pairwise_concordance", "numerator": "exact UTR matches", "denominator": "UTR comparisons", "truth_required": "yes for recall", "description": "Exact UTR-field agreement where UTR fields are available."},
        {"metric": "whole_model_concordance", "scope": "pairwise_concordance", "numerator": "intron-chain plus boundary matches", "denominator": "whole-model comparisons", "truth_required": "yes for recall", "description": "Combined structural model agreement; UTR-complete matches are labelled separately."},
    ])
    Path("benchmark_report.md").write_text(
        benchmark_report(models, whole_models, summaries, summary_rows, statuses),
        encoding="utf-8",
    )
    Path("benchmark_products_complete.txt").write_text(
        f"models={len(models)}\ncomparisons={len(comparisons)}\nwhole_models={len(whole_models)}\npairs={len(summaries)}\n",
        encoding="utf-8",
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest="command", required=True)
    pair = sub.add_parser("pair")
    pair.add_argument("--metadata", required=True)
    pair.add_argument("--source-a", required=True)
    pair.add_argument("--source-b", required=True)
    pair.add_argument("--evidence-a", default="")
    pair.add_argument("--evidence-b", default="")
    pair.add_argument("--intron", required=True)
    pair.add_argument("--ends", required=True)
    pair.add_argument("--utr", required=True)
    pair.add_argument("--structural", required=True)
    pair.add_argument("--output-models", required=True)
    pair.add_argument("--output-comparisons", required=True)
    pair.add_argument("--output-whole-models", required=True)
    pair.add_argument("--output-summary", required=True)
    pair.set_defaults(function=build_pair)

    combined = sub.add_parser("aggregate")
    combined.add_argument("--models", nargs="+", required=True)
    combined.add_argument("--comparisons", nargs="+", required=True)
    combined.add_argument("--whole-models", nargs="+", required=True)
    combined.add_argument("--summaries", nargs="+", required=True)
    combined.add_argument("--status", nargs="+", required=True)
    combined.set_defaults(function=aggregate)
    return root


if __name__ == "__main__":
    arguments = parser().parse_args()
    arguments.function(arguments)
