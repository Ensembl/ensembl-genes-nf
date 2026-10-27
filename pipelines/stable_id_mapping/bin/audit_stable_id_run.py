#!/usr/bin/env python3
"""Audit the main review buckets from a completed stable-ID mapper run."""

from __future__ import annotations

import argparse
import csv
import gzip
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from stable_id_mapping.gff3 import (
    GENE_FEATURE_TYPES,
    parent_ids,
    parse_attrs,
    split_stable_id,
)


BIOTYPE_KEYS = (
    "biotype",
    "gene_biotype",
    "transcript_biotype",
    "gbkey",
    "gene_type",
    "transcript_type",
)


@dataclass(frozen=True)
class GeneInfo:
    stable_id: str
    seqid: str
    start: str
    end: str
    strand: str
    biotype: str
    raw_id: str
    feature_type: str
    child_feature_types: tuple[str, ...] = ()

    @property
    def locus(self) -> str:
        return f"{self.seqid}:{self.start}-{self.end}({self.strand})"

    @property
    def annotation_class(self) -> str:
        if self.biotype:
            return self.biotype

        if self.feature_type and self.feature_type != "gene":
            return self.feature_type

        child_types = set(self.child_feature_types)
        if "mrna" in child_types:
            return "protein_coding_like_from_mRNA"
        if "transcript" in child_types:
            return "transcript_child_no_biotype"
        if child_types:
            return "+".join(sorted(child_types))
        return "<unknown>"


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarise and write review tables for missing genes, coordinate-only "
            "gene mappings, and new target genes."
        )
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Completed stable-ID mapper output directory",
    )
    parser.add_argument(
        "--decisions-tsv",
        type=Path,
        help="stable_id_decisions.tsv file to audit. Alternative to --run-dir.",
    )
    parser.add_argument(
        "--gene-locus-comparison",
        type=Path,
        help="lifton.gene_locus_comparison.tsv file to use.",
    )
    parser.add_argument("--ref-gff", type=Path, help="Reference GFF3 for biotypes/loci")
    parser.add_argument("--target-gff", type=Path, help="Target GFF3 for biotypes/loci")
    parser.add_argument(
        "--locus-comparison",
        type=Path,
        help=(
            "lifton.gene_locus_comparison.tsv to use. Defaults to "
            "RUN_DIR/matching/lifton.gene_locus_comparison.tsv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to RUN_DIR/reports/stable_id_audit",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Number of example rows to print per bucket",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()

    if args.run_dir is None and args.decisions_tsv is None:
        raise SystemExit("Provide either --run-dir or --decisions-tsv")

    output_dir = resolve_output_dir(args)
    summary = audit_run(
        run_dir=args.run_dir,
        ref_gff=args.ref_gff,
        target_gff=args.target_gff,
        output_dir=output_dir,
        locus_comparison=args.gene_locus_comparison or args.locus_comparison,
        decisions_tsv=args.decisions_tsv,
    )
    print_summary(summary, args.limit)


def resolve_output_dir(args: argparse.Namespace) -> Path:
    if args.output_dir is not None:
        return args.output_dir
    if args.run_dir is not None:
        return args.run_dir / "reports" / "stable_id_audit"
    return Path("stable_id_audit")


def audit_run(
    run_dir: Optional[Path],
    ref_gff: Optional[Path],
    target_gff: Optional[Path],
    output_dir: Path,
    locus_comparison: Optional[Path] = None,
    decisions_tsv: Optional[Path] = None,
) -> dict[str, object]:
    if decisions_tsv is not None:
        decisions_path = decisions_tsv
    elif run_dir is not None:
        decisions_path = run_dir / "stable_id_decisions.tsv"
    else:
        raise ValueError("Either run_dir or decisions_tsv is required")
    
    if not decisions_path.exists():
        raise FileNotFoundError(decisions_path)
    
    decisions = read_tsv(decisions_path)
    
    if locus_comparison is not None:
        locus_path = locus_comparison
    elif run_dir is not None:
        locus_path = run_dir / "matching" / "lifton.gene_locus_comparison.tsv"
    else:
        raise ValueError("Provide --gene-locus-comparison when using --decisions-tsv")
    locus_rows = read_tsv(locus_path)
    ref_genes = load_gene_info(ref_gff) if ref_gff else {}
    target_genes = load_gene_info(target_gff) if target_gff else {}

    locus_by_old = {
        core_id(row.get("lifton_gene", "")): row
        for row in locus_rows
        if core_id(row.get("lifton_gene", ""))
    }
    locus_usage = target_usage_from_locus_rows(locus_rows)

    gene_decisions = [row for row in decisions if row.get("type") == "gene"]
    gene_decisions_by_old = {
        row.get("old_stable_id", ""): row
        for row in gene_decisions
        if row.get("old_stable_id")
    }
    mapped_gene_decisions_by_target = {
        row.get("current_stable_id", ""): row
        for row in gene_decisions
        if row.get("action") == "mapped" and row.get("current_stable_id")
    }
    missing = [row for row in gene_decisions if row.get("action") == "missing"]
    coordinate_mapped = [
        row
        for row in gene_decisions
        if row.get("action") == "mapped" and "coordinate overlap" in row.get("reason", "")
    ]
    structural_mapped = [
        row
        for row in gene_decisions
        if row.get("action") == "mapped"
        and "lifton structural evidence" in row.get("reason", "")
    ]
    new = [row for row in gene_decisions if row.get("action") == "new"]

    output_dir.mkdir(parents=True, exist_ok=True)
    missing_rows = missing_gene_rows(
        missing,
        ref_genes,
        locus_by_old,
        mapped_gene_decisions_by_target,
    )
    coordinate_rows = coordinate_gene_rows(
        coordinate_mapped,
        ref_genes,
        target_genes,
        locus_by_old,
        mapped_gene_decisions_by_target,
    )
    new_rows = new_gene_rows(
        new,
        ref_genes,
        target_genes,
        locus_usage,
        gene_decisions_by_old,
    )

    paths = {
        "missing": output_dir / "missing_genes.tsv",
        "coordinate": output_dir / "coordinate_mapped_genes.tsv",
        "new": output_dir / "new_genes.tsv",
    }
    write_dicts(paths["missing"], missing_rows)
    write_dicts(paths["coordinate"], coordinate_rows)
    write_dicts(paths["new"], new_rows)

    return {
        "paths": paths,
        "locus_path": locus_path,
        "locus_rows_loaded": len(locus_rows),
        "feature_summary": build_feature_summary(decisions),
        "mapping_evidence_summary": build_mapping_evidence_summary(
            decisions
        ),
        "review_flags": {
            "coordinate_genes_below_0_50": sum(
                1
                for row in coordinate_mapped
                if to_float(row.get("score")) < 0.50
            ),
            "missing_genes_with_claimed_candidate": sum(
                1
                for row in missing_rows
                if (
                    row.get(
                        "target_gene_by_locus_claimed_by_old"
                    )
                    or row.get(
                        "structure_accepted_target_claimed_by_old"
                    )
                )
            ),
            "mapped_features_comparison_unavailable": sum(
                1
                for row in decisions
                if row.get("action") == "mapped"
                and "comparison unavailable"
                in row.get("reason", "")
            ),
        },
        "gene_annotation_summary": build_gene_annotation_summary(
            gene_decisions,
            ref_genes,
            target_genes,
        ),
        "counts": Counter(
            (row.get("type"), row.get("action"))
            for row in decisions
        ),
        "gene_counts": {
            "structural_mapped": len(structural_mapped),
            "coordinate_mapped": len(coordinate_mapped),
            "missing": len(missing),
            "new": len(new),
        },
        "missing_reason_counts": Counter(row.get("reason", "") for row in missing),
        "missing_locus_status": Counter(
            row.get("locus_vs_structure", "") or "<no locus row>"
            for row in missing_rows
        ),
        "coordinate_scores": tuple(to_float(row.get("score")) for row in coordinate_mapped),
        "coordinate_score_bands": Counter(
            score_band(to_float(row.get("score")))
            for row in coordinate_mapped
        ),
        "coordinate_locus_status": Counter(
            row.get("locus_vs_structure", "") or "<no locus row>"
            for row in coordinate_rows
        ),
        "new_annotation_classes": Counter(
            row.get("target_annotation_class", "<unknown>") for row in new_rows
        ),
        "new_target_id_in_ref": Counter(
            row.get("target_id_also_in_ref_gff", "no") for row in new_rows
        ),
        "new_target_id_in_ref_old_action": Counter(
            row.get("ref_same_id_action", "")
            for row in new_rows
            if row.get("target_id_also_in_ref_gff") == "yes"
        ),
        "new_with_locus_candidates": sum(
            1
            for row in new_rows
            if int(row.get("locus_candidate_for_old_count") or 0) > 0
        ),
        "new_with_structure_candidates": sum(
            1
            for row in new_rows
            if int(row.get("structure_candidate_for_old_count") or 0) > 0
        ),
        "missing_locus_target_claimed": Counter(
            yes_no(row.get("target_gene_by_locus_claimed_by_old"))
            for row in missing_rows
        ),
        "missing_locus_status_and_claimed": Counter(
            (
                row.get("locus_vs_structure", "") or "<no locus row>",
                yes_no(row.get("target_gene_by_locus_claimed_by_old")),
            )
            for row in missing_rows
        ),
        "missing_structure_target_claimed": Counter(
            yes_no(row.get("structure_accepted_target_claimed_by_old"))
            for row in missing_rows
        ),
        "missing_rows": missing_rows,
        "coordinate_rows": coordinate_rows,
        "new_rows": new_rows,
    }


def missing_gene_rows(
    rows: list[dict[str, str]],
    ref_genes: dict[str, GeneInfo],
    locus_by_old: dict[str, dict[str, str]],
    mapped_gene_decisions_by_target: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        old_id = row.get("old_stable_id", "")
        old = ref_genes.get(old_id)
        locus = locus_by_old.get(old_id, {})
        target_gene_by_locus = core_id(locus.get("target_gene_by_locus", ""))
        structural_target_gene = core_id(locus.get("structure_accepted_target_gene", ""))
        best_tx_target_gene = core_id(locus.get("best_tx_structure_target_gene", ""))
        locus_target_claim = mapped_gene_decisions_by_target.get(target_gene_by_locus, {})
        structural_target_claim = mapped_gene_decisions_by_target.get(
            structural_target_gene,
            {},
        )
        best_tx_target_claim = mapped_gene_decisions_by_target.get(best_tx_target_gene, {})
        out.append(
            {
                "old_stable_id": old_id,
                "old_biotype": old.biotype if old else "",
                "old_locus": old.locus if old else "",
                "reason": row.get("reason", ""),
                "locus_vs_structure": locus.get("locus_vs_structure", ""),
                "locus_score": locus.get("locus_score", ""),
                "target_gene_by_locus": target_gene_by_locus,
                "target_gene_by_locus_claimed_by_old": locus_target_claim.get(
                    "old_stable_id",
                    "",
                ),
                "target_gene_by_locus_claim_score": locus_target_claim.get("score", ""),
                "structure_accepted_target_gene": structural_target_gene,
                "structure_accepted_target_claimed_by_old": structural_target_claim.get(
                    "old_stable_id",
                    "",
                ),
                "structure_accepted_target_claim_score": structural_target_claim.get(
                    "score",
                    "",
                ),
                "best_tx_structure_target_gene": best_tx_target_gene,
                "best_tx_structure_target_claimed_by_old": best_tx_target_claim.get(
                    "old_stable_id",
                    "",
                ),
                "best_tx_structure_target_claim_score": best_tx_target_claim.get(
                    "score",
                    "",
                ),
                "best_tx_structure_score": locus.get("best_tx_structure_score", ""),
                "lifton_locus": format_locus(
                    locus.get("lifton_contig", ""),
                    locus.get("lifton_start", ""),
                    locus.get("lifton_end", ""),
                    locus.get("lifton_strand", ""),
                ),
            }
        )
    return sorted(
        out,
        key=lambda item: (
            -to_float(item.get("locus_score")),
            item.get("reason", ""),
            item.get("old_stable_id", ""),
        ),
    )


def coordinate_gene_rows(
    rows: list[dict[str, str]],
    ref_genes: dict[str, GeneInfo],
    target_genes: dict[str, GeneInfo],
    locus_by_old: dict[str, dict[str, str]],
    mapped_gene_decisions_by_target: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        old_id = row.get("old_stable_id", "")
        target_id = row.get("current_stable_id", "")
        old = ref_genes.get(old_id)
        target = target_genes.get(target_id)
        locus = locus_by_old.get(old_id, {})
        target_gene_by_locus = core_id(locus.get("target_gene_by_locus", ""))
        structural_target_gene = core_id(locus.get("structure_accepted_target_gene", ""))
        best_tx_target_gene = core_id(locus.get("best_tx_structure_target_gene", ""))
        structural_target_claim = mapped_gene_decisions_by_target.get(
            structural_target_gene,
            {},
        )
        best_tx_target_claim = mapped_gene_decisions_by_target.get(best_tx_target_gene, {})
        out.append(
            {
                "old_stable_id": old_id,
                "target_stable_id": target_id,
                "score": row.get("score", ""),
                "old_biotype": old.biotype if old else "",
                "old_annotation_class": old.annotation_class if old else "",
                "target_biotype": target.biotype if target else "",
                "target_annotation_class": target.annotation_class if target else "",
                "old_locus": old.locus if old else "",
                "target_locus": target.locus if target else "",
                "reason": row.get("reason", ""),
                "locus_vs_structure": locus.get("locus_vs_structure", ""),
                "locus_score": locus.get("locus_score", ""),
                "target_gene_by_locus": target_gene_by_locus,
                "structure_accepted_target_gene": structural_target_gene,
                "structure_accepted_target_claimed_by_old": structural_target_claim.get(
                    "old_stable_id",
                    "",
                ),
                "structure_accepted_target_claim_score": structural_target_claim.get(
                    "score",
                    "",
                ),
                "best_tx_structure_target_gene": best_tx_target_gene,
                "best_tx_structure_target_claimed_by_old": best_tx_target_claim.get(
                    "old_stable_id",
                    "",
                ),
                "best_tx_structure_target_claim_score": best_tx_target_claim.get(
                    "score",
                    "",
                ),
                "best_tx_structure_score": locus.get("best_tx_structure_score", ""),
            }
        )
    return sorted(
        out,
        key=lambda item: (
            to_float(item.get("score")),
            -to_float(item.get("locus_score")),
            item.get("old_stable_id", ""),
        ),
    )


def new_gene_rows(
    rows: list[dict[str, str]],
    ref_genes: dict[str, GeneInfo],
    target_genes: dict[str, GeneInfo],
    locus_usage: dict[str, dict[str, object]],
    gene_decisions_by_old: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        target_id = row.get("current_stable_id", "")
        target = target_genes.get(target_id)
        usage = locus_usage.get(target_id, {})
        same_id_old_decision = gene_decisions_by_old.get(target_id, {})
        out.append(
            {
                "target_stable_id": target_id,
                "assigned_new_stable_id": row.get("new_stable_id", ""),
                "target_biotype": target.biotype if target else "",
                "target_annotation_class": target.annotation_class if target else "",
                "target_locus": target.locus if target else "",
                "target_id_also_in_ref_gff": "yes" if target_id in ref_genes else "no",
                "ref_same_id_action": same_id_old_decision.get("action", ""),
                "ref_same_id_mapped_target": same_id_old_decision.get("current_stable_id", ""),
                "ref_same_id_score": same_id_old_decision.get("score", ""),
                "ref_same_id_reason": same_id_old_decision.get("reason", ""),
                "locus_candidate_for_old_count": str(
                    usage.get("locus_candidate_count", 0)
                ),
                "best_old_by_locus": str(usage.get("best_old_by_locus", "")),
                "best_locus_score": str(usage.get("best_locus_score", "")),
                "structure_candidate_for_old_count": str(
                    usage.get("structure_candidate_count", 0)
                ),
                "best_old_by_structure": str(usage.get("best_old_by_structure", "")),
                "best_structure_score": str(usage.get("best_structure_score", "")),
                "reason": row.get("reason", ""),
            }
        )
    return sorted(
        out,
        key=lambda item: (
            item.get("target_biotype", ""),
            -int(item.get("locus_candidate_for_old_count") or 0),
            item.get("target_stable_id", ""),
        ),
    )


def target_usage_from_locus_rows(
    locus_rows: list[dict[str, str]],
) -> dict[str, dict[str, object]]:
    usage: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "locus_candidate_count": 0,
            "best_locus_score": "",
            "best_old_by_locus": "",
            "structure_candidate_count": 0,
            "best_structure_score": "",
            "best_old_by_structure": "",
        }
    )
    for row in locus_rows:
        old_id = core_id(row.get("lifton_gene", ""))
        locus_target = core_id(row.get("target_gene_by_locus", ""))
        if locus_target:
            item = usage[locus_target]
            item["locus_candidate_count"] = int(item["locus_candidate_count"]) + 1
            locus_score = to_float(row.get("locus_score"))
            if locus_score > to_float(str(item["best_locus_score"])):
                item["best_locus_score"] = row.get("locus_score", "")
                item["best_old_by_locus"] = old_id

        structure_target = core_id(row.get("structure_accepted_target_gene", ""))
        if structure_target:
            item = usage[structure_target]
            item["structure_candidate_count"] = (
                int(item["structure_candidate_count"]) + 1
            )
            structure_score = to_float(row.get("structure_weighted_score"))
            if structure_score > to_float(str(item["best_structure_score"])):
                item["best_structure_score"] = row.get("structure_weighted_score", "")
                item["best_old_by_structure"] = old_id
    return dict(usage)


def load_gene_info(path: Path) -> dict[str, GeneInfo]:
    gene_rows: dict[str, tuple[str, str, str, str, str, str, str]] = {}
    child_types_by_gene: dict[str, set[str]] = defaultdict(set)
    for fields in iter_gff3_rows(path):
        seqid, _source, feature_type, start, end, _score, strand, _phase, attrs_text = fields
        feature_type_lc = feature_type.lower()
        attrs = parse_attrs(attrs_text)
        if feature_type_lc in GENE_FEATURE_TYPES:
            stable_id, _version = split_stable_id(attrs.get("ID"))
            if not stable_id:
                continue
            gene_rows[stable_id] = (
                seqid,
                start,
                end,
                strand,
                first_attr(attrs, BIOTYPE_KEYS) or "",
                attrs.get("ID", ""),
                feature_type_lc,
            )
            continue
        for parent_id, _version in parent_ids(attrs.get("Parent")):
            child_types_by_gene[parent_id].add(feature_type_lc)

    return {
        stable_id: GeneInfo(
            stable_id=stable_id,
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            biotype=biotype,
            raw_id=raw_id,
            feature_type=feature_type,
            child_feature_types=tuple(sorted(child_types_by_gene.get(stable_id, ()))),
        )
        for stable_id, (
            seqid,
            start,
            end,
            strand,
            biotype,
            raw_id,
            feature_type,
        ) in gene_rows.items()
    }


def iter_gff3_rows(path: Path) -> Iterable[list[str]]:
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9:
                yield fields


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_dicts(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        fieldnames = list(rows[0])
    else:
        fieldnames = ["note"]
        rows = [{"note": "no rows"}]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def percentage(part: int, total: int) -> Optional[float]:
    if total == 0:
        return None
    return (part / total) * 100.0


def build_feature_summary(
    decisions: list[dict[str, str]],
) -> list[dict[str, object]]:
    feature_order = (
        "gene",
        "transcript",
        "translation",
        "exon",
    )
    counts = Counter(
        (row.get("type", ""), row.get("action", ""))
        for row in decisions
    )

    rows: list[dict[str, object]] = []

    for feature_type in feature_order:
        retained = counts[(feature_type, "mapped")]
        missing = counts[(feature_type, "missing")]
        new = counts[(feature_type, "new")]

        reference_total = retained + missing
        target_total = retained + new

        mapped_rows = [
            row
            for row in decisions
            if row.get("type") == feature_type
            and row.get("action") == "mapped"
        ]

        version_unchanged = sum(
            1
            for row in mapped_rows
            if int(row.get("new_version") or 0)
            == int(row.get("old_version") or 0)
        )
        version_incremented = sum(
            1
            for row in mapped_rows
            if int(row.get("new_version") or 0)
            > int(row.get("old_version") or 0)
        )
        comparison_unavailable = sum(
            1
            for row in mapped_rows
            if "comparison unavailable" in row.get("reason", "")
        )

        rows.append(
            {
                "feature_type": feature_type,
                "reference_total": reference_total,
                "retained": retained,
                "retained_percent": percentage(
                    retained,
                    reference_total,
                ),
                "missing": missing,
                "target_total": target_total,
                "new": new,
                "new_percent": percentage(
                    new,
                    target_total,
                ),
                "version_unchanged": version_unchanged,
                "version_incremented": version_incremented,
                "comparison_unavailable": comparison_unavailable,
            }
        )

    return rows


def build_mapping_evidence_summary(
    decisions: list[dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for feature_type in ("gene", "transcript"):
        mapped_rows = [
            row
            for row in decisions
            if row.get("type") == feature_type
            and row.get("action") == "mapped"
        ]

        total = len(mapped_rows)
        structural = sum(
            1
            for row in mapped_rows
            if "lifton structural evidence"
            in row.get("reason", "")
        )
        coordinate = sum(
            1
            for row in mapped_rows
            if "coordinate overlap"
            in row.get("reason", "")
        )
        other = total - structural - coordinate

        rows.append(
            {
                "feature_type": feature_type,
                "total": total,
                "structural": structural,
                "structural_percent": percentage(
                    structural,
                    total,
                ),
                "coordinate": coordinate,
                "coordinate_percent": percentage(
                    coordinate,
                    total,
                ),
                "other": other,
                "other_percent": percentage(
                    other,
                    total,
                ),
            }
        )

    return rows


def gene_annotation_key(
    gene_info: Optional[GeneInfo],
) -> tuple[str, str]:
    if gene_info is None:
        return "<unknown>", "<unknown>"

    return (
        gene_info.feature_type,
        gene_info.annotation_class,
    )


def build_gene_annotation_summary(
    gene_decisions: list[dict[str, str]],
    ref_genes: dict[str, GeneInfo],
    target_genes: dict[str, GeneInfo],
) -> dict[str, object]:
    reference_counts: dict[
        tuple[str, str],
        Counter,
    ] = defaultdict(Counter)

    target_counts: dict[
        tuple[str, str],
        Counter,
    ] = defaultdict(Counter)

    transitions: Counter = Counter()

    for row in gene_decisions:
        action = row.get("action", "")

        old_info = ref_genes.get(
            row.get("old_stable_id", "")
        )
        target_info = target_genes.get(
            row.get("current_stable_id", "")
        )

        if action in {"mapped", "missing"}:
            reference_key = gene_annotation_key(old_info)
            reference_counts[reference_key]["total"] += 1
            reference_counts[reference_key][action] += 1

        if action in {"mapped", "new"}:
            target_key = gene_annotation_key(target_info)
            target_counts[target_key]["total"] += 1
            target_counts[target_key][action] += 1

        if action == "mapped":
            reference_key = gene_annotation_key(old_info)
            target_key = gene_annotation_key(target_info)

            if reference_key != target_key:
                transitions[
                    (
                        reference_key[0],
                        reference_key[1],
                        target_key[0],
                        target_key[1],
                    )
                ] += 1

    reference_rows: list[dict[str, object]] = []

    for (
        feature_type,
        annotation_class,
    ), values in sorted(
        reference_counts.items(),
        key=lambda item: (
            -item[1]["total"],
            item[0][0],
            item[0][1],
        ),
    ):
        total = values["total"]
        retained = values["mapped"]
        missing = values["missing"]

        reference_rows.append(
            {
                "feature_type": feature_type,
                "annotation_class": annotation_class,
                "reference_total": total,
                "retained": retained,
                "retained_percent": percentage(
                    retained,
                    total,
                ),
                "missing": missing,
            }
        )

    target_rows: list[dict[str, object]] = []

    for (
        feature_type,
        annotation_class,
    ), values in sorted(
        target_counts.items(),
        key=lambda item: (
            -item[1]["total"],
            item[0][0],
            item[0][1],
        ),
    ):
        total = values["total"]
        retained = values["mapped"]
        new = values["new"]

        target_rows.append(
            {
                "feature_type": feature_type,
                "annotation_class": annotation_class,
                "target_total": total,
                "retained": retained,
                "new": new,
                "new_percent": percentage(
                    new,
                    total,
                ),
            }
        )

    transition_rows = [
        {
            "reference_feature_type": reference_feature_type,
            "reference_annotation_class": reference_annotation_class,
            "target_feature_type": target_feature_type,
            "target_annotation_class": target_annotation_class,
            "count": count,
        }
        for (
            reference_feature_type,
            reference_annotation_class,
            target_feature_type,
            target_annotation_class,
        ), count in transitions.most_common()
    ]

    return {
        "reference": reference_rows,
        "target": target_rows,
        "transitions": transition_rows,
    }


def format_percentage(value: object) -> str:
    if value is None:
        return "-"
    return f"{float(value):.2f}%"


def format_count(
    value: object,
    applicable: bool = True,
) -> str:
    if not applicable:
        return "-"
    return f"{int(value):,}"


def print_feature_summary(
    rows: list[dict[str, object]],
) -> None:
    print("STABLE-ID MAPPING SUMMARY")
    print()
    print(
        f"{'Feature':<12}"
        f"{'Reference':>12}"
        f"{'Retained':>12}"
        f"{'Retained %':>13}"
        f"{'Missing':>12}"
        f"{'Target':>12}"
        f"{'New':>12}"
        f"{'New %':>11}"
    )

    for row in rows:
        reference_total = int(row["reference_total"])
        target_total = int(row["target_total"])
        reference_applicable = reference_total > 0
        target_applicable = target_total > 0

        print(
            f"{str(row['feature_type']):<12}"
            f"{format_count(reference_total, reference_applicable):>12}"
            f"{format_count(row['retained'], reference_applicable):>12}"
            f"{format_percentage(row['retained_percent']):>13}"
            f"{format_count(row['missing'], reference_applicable):>12}"
            f"{format_count(target_total, target_applicable):>12}"
            f"{format_count(row['new'], target_applicable):>12}"
            f"{format_percentage(row['new_percent']):>11}"
        )

    print()
    print("VERSION OUTCOMES AMONG RETAINED IDS")
    print()
    print(
        f"{'Feature':<12}"
        f"{'Retained':>12}"
        f"{'Unchanged':>12}"
        f"{'Unchanged %':>13}"
        f"{'Incremented':>13}"
        f"{'Incremented %':>15}"
        f"{'Unavailable':>13}"
    )

    for row in rows:
        retained = int(row["retained"])

        if retained == 0:
            continue

        unchanged = int(row["version_unchanged"])
        incremented = int(row["version_incremented"])

        print(
            f"{str(row['feature_type']):<12}"
            f"{retained:>12,}"
            f"{unchanged:>12,}"
            f"{format_percentage(percentage(unchanged, retained)):>13}"
            f"{incremented:>13,}"
            f"{format_percentage(percentage(incremented, retained)):>15}"
            f"{int(row['comparison_unavailable']):>13,}"
        )

    print()
    print(
        "Retained % = mapped reference features / all reference features."
    )
    print(
        "New % = newly assigned target features / all target features."
    )
    print(
        "Exons are target-only assignments and therefore have no "
        "retention percentage."
    )


def print_gene_annotation_summary(
    summary: dict[str, object],
    transition_limit: int,
) -> None:
    reference_rows = summary["reference"]
    target_rows = summary["target"]
    transitions = summary["transitions"]

    assert isinstance(reference_rows, list)
    assert isinstance(target_rows, list)
    assert isinstance(transitions, list)

    print()
    print("REFERENCE GENE CLASSES")
    print()
    print(
        f"{'GFF type / annotation class':<45}"
        f"{'Reference':>12}"
        f"{'Retained':>12}"
        f"{'Retained %':>13}"
        f"{'Missing':>12}"
    )

    for row in reference_rows:
        label = (
            f"{row['feature_type']} / "
            f"{row['annotation_class']}"
        )
        print(
            f"{label:<45}"
            f"{int(row['reference_total']):>12,}"
            f"{int(row['retained']):>12,}"
            f"{format_percentage(row['retained_percent']):>13}"
            f"{int(row['missing']):>12,}"
        )

    print()
    print("TARGET GENE CLASSES")
    print()
    print(
        f"{'GFF type / annotation class':<45}"
        f"{'Target':>12}"
        f"{'Retained':>12}"
        f"{'New':>12}"
        f"{'New %':>11}"
    )

    for row in target_rows:
        label = (
            f"{row['feature_type']} / "
            f"{row['annotation_class']}"
        )
        print(
            f"{label:<45}"
            f"{int(row['target_total']):>12,}"
            f"{int(row['retained']):>12,}"
            f"{int(row['new']):>12,}"
            f"{format_percentage(row['new_percent']):>11}"
        )

    print()
    print("ANNOTATION-CLASS CHANGES AMONG RETAINED GENES")
    print()

    if not transitions:
        print("  none")
        return

    total_retained = sum(
        int(row["retained"])
        for row in reference_rows
    )

    for row in transitions[:transition_limit]:
        reference_label = (
            f"{row['reference_feature_type']} / "
            f"{row['reference_annotation_class']}"
        )
        target_label = (
            f"{row['target_feature_type']} / "
            f"{row['target_annotation_class']}"
        )
        count = int(row["count"])

        print(
            f"  {reference_label} -> {target_label}: "
            f"{count:,} "
            f"({format_percentage(percentage(count, total_retained))})"
        )


def print_mapping_evidence_summary(
    rows: list[dict[str, object]],
) -> None:
    print()
    print("MAPPING EVIDENCE AMONG RETAINED IDS")
    print()
    print(
        f"{'Feature':<12}"
        f"{'Mapped':>12}"
        f"{'Structural':>12}"
        f"{'Structural %':>14}"
        f"{'Coordinate':>12}"
        f"{'Coordinate %':>14}"
        f"{'Other':>10}"
        f"{'Other %':>10}"
    )

    for row in rows:
        print(
            f"{str(row['feature_type']):<12}"
            f"{int(row['total']):>12,}"
            f"{int(row['structural']):>12,}"
            f"{format_percentage(row['structural_percent']):>14}"
            f"{int(row['coordinate']):>12,}"
            f"{format_percentage(row['coordinate_percent']):>14}"
            f"{int(row['other']):>10,}"
            f"{format_percentage(row['other_percent']):>10}"
        )


def print_review_flags(
    flags: dict[str, object],
) -> None:
    print()
    print("REVIEW FLAGS")
    print()
    print(
        "  Coordinate-only genes with score below 0.50: "
        f"{int(flags['coordinate_genes_below_0_50']):,}"
    )
    print(
        "  Missing genes with a candidate claimed elsewhere: "
        f"{int(flags['missing_genes_with_claimed_candidate']):,}"
    )
    print(
        "  Mapped features with unavailable version comparisons: "
        f"{int(flags['mapped_features_comparison_unavailable']):,}"
    )


def print_summary(summary: dict[str, object], limit: int) -> None:
    feature_summary = summary["feature_summary"]
    assert isinstance(feature_summary, list)
    print_feature_summary(feature_summary)

    mapping_evidence_summary = summary[
        "mapping_evidence_summary"
    ]
    assert isinstance(mapping_evidence_summary, list)
    print_mapping_evidence_summary(
        mapping_evidence_summary
    )

    review_flags = summary["review_flags"]
    assert isinstance(review_flags, dict)
    print_review_flags(review_flags)

    gene_annotation_summary = summary["gene_annotation_summary"]
    assert isinstance(gene_annotation_summary, dict)
    print_gene_annotation_summary(
        gene_annotation_summary,
        transition_limit=limit,
    )

    print()
    print("DETAILED MAPPING AUDIT")
    print()

    decision_counts = summary["counts"]
    assert isinstance(decision_counts, Counter)

    print("All stable-ID decisions:")
    for (feature_type, action), count in sorted(decision_counts.items()):
        print(f"  {feature_type}:{action}: {count}")
    
    gene_counts = summary["gene_counts"]
    assert isinstance(gene_counts, dict)
    print("Gene decision audit:")
    for key in ("structural_mapped", "coordinate_mapped", "missing", "new"):
        print(f"  {key}: {gene_counts[key]}")
    print(
        "Gene locus comparison rows loaded: "
        f"{summary['locus_rows_loaded']} ({summary['locus_path']})"
    )

    coordinate_scores = summary["coordinate_scores"]
    assert isinstance(coordinate_scores, tuple)
    if coordinate_scores:
        print(
            "Coordinate-only gene scores: "
            f"min={min(coordinate_scores):.3f}, "
            f"median={statistics.median(coordinate_scores):.3f}, "
            f"max={max(coordinate_scores):.3f}"
        )
    print_counter("Coordinate-only score bands", summary["coordinate_score_bands"])

    print_counter("Missing gene reasons", summary["missing_reason_counts"])
    print_counter("Missing gene locus status", summary["missing_locus_status"])
    print_counter(
        "Missing gene locus candidate target already claimed",
        summary["missing_locus_target_claimed"],
    )
    print_counter(
        "Missing gene structural target already claimed",
        summary["missing_structure_target_claimed"],
    )
    print_pair_counter(
        "Missing gene locus status + target claimed",
        summary["missing_locus_status_and_claimed"],
    )
    print_counter("Coordinate-only gene locus status", summary["coordinate_locus_status"])
    print_counter(
        "New gene annotation classes",
        summary["new_annotation_classes"],
        limit=10,
    )
    print_counter("New target ID also present in ref GFF", summary["new_target_id_in_ref"])
    print_counter(
        "New target IDs also in ref GFF: old-ID decision action",
        summary["new_target_id_in_ref_old_action"],
    )
    print(
        "New genes seen as candidate targets: "
        f"locus={summary['new_with_locus_candidates']}, "
        f"structure={summary['new_with_structure_candidates']}"
    )

    paths = summary["paths"]
    assert isinstance(paths, dict)
    print("Audit tables:")
    for key in ("missing", "coordinate", "new"):
        print(f"  {key}: {paths[key]}")

    print_examples("Missing genes", summary["missing_rows"], limit)
    print_examples("Coordinate-only mapped genes", summary["coordinate_rows"], limit)
    print_examples("New target genes", summary["new_rows"], limit)


def print_counter(label: str, value: object, limit: int = 5) -> None:
    assert isinstance(value, Counter)
    if not value:
        print(f"{label}: none")
        return
    print(f"{label}:")
    for key, count in value.most_common(limit):
        print(f"  {key or '<blank>'}: {count}")


def print_pair_counter(label: str, value: object, limit: int = 10) -> None:
    assert isinstance(value, Counter)
    if not value:
        print(f"{label}: none")
        return
    print(f"{label}:")
    for (left, right), count in value.most_common(limit):
        print(f"  {left}, claimed={right}: {count}")


def print_examples(label: str, value: object, limit: int) -> None:
    assert isinstance(value, list)
    print(f"{label} examples:")
    if not value:
        print("  none")
        return
    for row in value[:limit]:
        if "assigned_new_stable_id" in row:
            visible = [
                row.get("target_stable_id", ""),
                row.get("assigned_new_stable_id", ""),
                row.get("target_annotation_class", ""),
                (
                    "in_ref="
                    f"{row.get('target_id_also_in_ref_gff', '')};"
                    f"old_action={row.get('ref_same_id_action', '')}"
                ),
                row.get("reason", ""),
            ]
        elif "target_stable_id" in row:
            visible = [
                row.get("old_stable_id", ""),
                row.get("target_stable_id", ""),
                row.get("score", ""),
                row.get("locus_vs_structure", ""),
                row.get("reason", ""),
            ]
        else:
            visible = [
                row.get("old_stable_id", ""),
                row.get("old_biotype", ""),
                row.get("locus_score", ""),
                row.get("locus_vs_structure", ""),
                (
                    "locus_claimed_by="
                    f"{row.get('target_gene_by_locus_claimed_by_old', '')};"
                    "structure_claimed_by="
                    f"{row.get('structure_accepted_target_claimed_by_old', '')}"
                ),
            ]
        print("  " + "\t".join(str(item) for item in visible))


def core_id(value: Optional[str]) -> str:
    stable_id, _version = split_stable_id(value)
    return stable_id or ""


def first_attr(attrs: dict[str, str], keys: Iterable[str]) -> Optional[str]:
    for key in keys:
        value = attrs.get(key)
        if value:
            return value
    return None


def format_locus(seqid: str, start: str, end: str, strand: str) -> str:
    if not seqid:
        return ""
    return f"{seqid}:{start}-{end}({strand})"


def to_float(value: Optional[str]) -> float:
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def yes_no(value: Optional[str]) -> str:
    return "yes" if value else "no"


def score_band(score: float) -> str:
    if score < 0.5:
        return "<0.50"
    if score < 0.75:
        return "0.50-0.75"
    if score < 0.95:
        return "0.75-0.95"
    return ">=0.95"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
