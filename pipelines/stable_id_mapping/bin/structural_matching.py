#!/usr/bin/env python3
"""Run LiftOn structural matching for one species."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from stable_id_mapping.lifton_matching import LiftonMatchConfig, run_lifton_matching
from stable_id_mapping.rules import DEFAULT_RULES_PATH, load_mapping_rules


@dataclass(frozen=True)
class StructuralMatchingConfig:
    lifton_gff: Path
    target_gff: Path
    out_prefix: Path
    rules_config: Path
    output_json: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lifton-gff", type=Path, required=True)
    parser.add_argument("--target-gff", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--rules-config", type=Path, default=DEFAULT_RULES_PATH)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    rules = load_mapping_rules(args.rules_config).structural_matching

    summary = run_lifton_matching(
        LiftonMatchConfig(
            lifton_gff=args.lifton_gff,
            target_gff=args.target_gff,
            out_prefix=args.out_prefix,
            window=rules.window,
            topk=rules.topk,
            min_score=rules.min_score,
            good=rules.good_score,
            confident=rules.confident_score,
            gene_fraction=rules.gene_fraction,
            score_weights=dict(rules.score_weights),
        )
    )

    args.output_json.write_text(
        json.dumps(
            {
                "transcript_pairs": summary.transcript_pairs,
                "gene_pairs": summary.gene_pairs,
                "transcript_pairs_path": str(summary.transcript_pairs_path),
                "gene_pairs_path": str(summary.gene_pairs_path),
                "gene_locus_comparison_path": str(summary.gene_locus_comparison_path),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
