#!/usr/bin/env python3
"""Write final one-to-one relationships between two annotation GFF3 files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from stable_id_mapping.decisions import build_decisions
from stable_id_mapping.gff3 import parse_gff3, parse_missing_gene_ids
from stable_id_mapping.ids import StableIdRange, collect_reserved_ids, make_allocator
from stable_id_mapping.rules import DEFAULT_RULES_PATH, load_mapping_rules
from stable_id_mapping.scoring import load_lifton_score_evidence


COMPARISON_RANGES = {
    "gene": StableIdRange("COMPAREG", 1, 999_999_999_999, 12),
    "transcript": StableIdRange("COMPARET", 1, 999_999_999_999, 12),
    "translation": StableIdRange("COMPAREP", 1, 999_999_999_999, 12),
    "exon": StableIdRange("COMPAREE", 1, 999_999_999_999, 12),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument("--target-gff", type=Path, required=True)
    parser.add_argument("--mapped-gff", type=Path, required=True)
    parser.add_argument("--missing-report", type=Path, required=True)
    parser.add_argument("--transcript-pairs", type=Path, required=True)
    parser.add_argument("--gene-pairs", type=Path, required=True)
    parser.add_argument("--rules-config", type=Path, default=DEFAULT_RULES_PATH)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--no-translations", action="store_true")
    return parser.parse_args()


def relationship_method(reason: str) -> str:
    if "lifton structural evidence" in reason:
        return "structural"
    if "coordinate overlap" in reason:
        return "coordinate"
    if "single old and target translation" in reason:
        return "single_translation_pair"
    if "parent_transcript_lifton_score" in reason:
        return "parent_transcript"
    return "derived"


def run_comparison(
    ref_gff: Path,
    target_gff: Path,
    mapped_gff: Path,
    missing_report: Path,
    transcript_pairs: Path,
    gene_pairs: Path,
    rules_config: Path,
    output_tsv: Path,
    include_translations: bool,
) -> None:
    rules = load_mapping_rules(rules_config)
    ref_features = parse_gff3(ref_gff)
    target_features = parse_gff3(target_gff, require_exon_ids=True)
    mapped_features = parse_gff3(mapped_gff)
    score_evidence = load_lifton_score_evidence(transcript_pairs, gene_pairs)
    reserved_ids = collect_reserved_ids(
        ref_features,
        target_features,
        mapped_features,
    )
    allocators = {
        feature_type: make_allocator(id_range, reserved_ids)
        for feature_type, id_range in COMPARISON_RANGES.items()
    }

    decisions = build_decisions(
        ref_features,
        target_features,
        mapped_features,
        parse_missing_gene_ids(missing_report),
        allocators,
        mapping_session_id=0,
        min_overlap=rules.coordinate_overlap.min_overlap,
        include_translations=include_translations,
        score_evidence=score_evidence,
    )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "feature_type",
                "reference_stable_id",
                "target_stable_id",
                "method",
                "score",
                "reason",
            ]
        )
        for decision in decisions:
            if decision.action != "mapped":
                continue
            writer.writerow(
                [
                    decision.feature_type,
                    decision.old_stable_id,
                    decision.current_stable_id,
                    relationship_method(decision.reason),
                    f"{decision.score:.6f}",
                    decision.reason,
                ]
            )


def main() -> None:
    args = parse_args()
    run_comparison(
        ref_gff=args.ref_gff,
        target_gff=args.target_gff,
        mapped_gff=args.mapped_gff,
        missing_report=args.missing_report,
        transcript_pairs=args.transcript_pairs,
        gene_pairs=args.gene_pairs,
        rules_config=args.rules_config,
        output_tsv=args.output_tsv,
        include_translations=not args.no_translations,
    )


if __name__ == "__main__":
    main()
