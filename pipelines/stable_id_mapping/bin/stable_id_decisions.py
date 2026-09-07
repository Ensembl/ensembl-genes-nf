#!/usr/bin/env python3
"""Generate stable-ID mapping decisions for one species."""

from __future__ import annotations

import argparse
import json
import tempfile
from collections import Counter
from pathlib import Path

from stable_id_mapping.config import StableIdEventConfig
from stable_id_mapping.ids import parse_id_range
from stable_id_mapping.pipeline import run_pipeline
from stable_id_mapping.rules import DEFAULT_RULES_PATH, load_mapping_rules
from stable_id_mapping.scoring import load_lifton_score_evidence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument("--target-gff", type=Path, required=True)
    parser.add_argument("--ref-fasta", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--mapped-gff", type=Path, required=True)
    parser.add_argument("--missing-report", type=Path, required=True)
    parser.add_argument("--transcript-pairs", type=Path, required=True)
    parser.add_argument("--gene-pairs", type=Path, required=True)
    parser.add_argument("--mapping-session-id", type=int, required=True)
    parser.add_argument("--gene-range", type=parse_id_range, required=True)
    parser.add_argument("--transcript-range", type=parse_id_range, required=True)
    parser.add_argument("--translation-range", type=parse_id_range, required=True)
    parser.add_argument("--decisions-tsv", type=Path, required=True)
    parser.add_argument("--score-evidence-tsv", type=Path, required=True)
    parser.add_argument("--rules-config", type=Path, default=DEFAULT_RULES_PATH)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--no-translations", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rules = load_mapping_rules(args.rules_config)
    score_evidence = load_lifton_score_evidence(
        args.transcript_pairs,
        args.gene_pairs,
    )

    # SQL generation remains a later Nextflow stage. The current pipeline API
    # requires an SQL path, so discard that intermediate artifact here.
    with tempfile.TemporaryDirectory(prefix="stable_id_decisions_") as temp_dir:
        decisions = run_pipeline(
            StableIdEventConfig(
                ref_gff=args.ref_gff,
                target_gff=args.target_gff,
                ref_fasta=args.ref_fasta,
                target_fasta=args.target_fasta,
                mapped_gff=args.mapped_gff,
                report=args.missing_report,
                mapping_session_id=args.mapping_session_id,
                gene_range=args.gene_range,
                transcript_range=args.transcript_range,
                translation_range=args.translation_range,
                output_sql=Path(temp_dir) / "deferred.sql",
                output_tsv=args.decisions_tsv,
                db_name=args.db_name,
                include_translations=not args.no_translations,
                min_overlap=rules.coordinate_overlap.min_overlap,
                score_evidence=score_evidence,
                output_score_evidence_tsv=args.score_evidence_tsv,
            )
        )

    action_counts = Counter(decision.action for decision in decisions)
    feature_counts = Counter(decision.feature_type for decision in decisions)
    args.output_json.write_text(
        json.dumps(
            {
                "db_name": args.db_name,
                "decisions": len(decisions),
                "assigned": sum(
                    decision.new_stable_id is not None for decision in decisions
                ),
                "missing": action_counts["missing"],
                "exons_assigned": feature_counts["exon"],
                "actions": dict(sorted(action_counts.items())),
                "feature_types": dict(sorted(feature_counts.items())),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
