#!/usr/bin/env python3
"""Run LiftOn projection for one species and write the missing-gene report.

Thin wrapper around stable_id_mapping.lifton.run_lifton and
stable_id_mapping.reports.write_missing_gene_report — reuses the existing
pipeline package rather than reimplementing any matching/decision logic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from stable_id_mapping.lifton import DEFAULT_LIFTON_FEATURE_TYPES, LiftonRunConfig, run_lifton
from stable_id_mapping.reports import write_missing_gene_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-fasta", type=Path, required=True)
    parser.add_argument("--ref-gff", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--output-gff", type=Path, required=True)
    parser.add_argument("--missing-report", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--executable", default="lifton")
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument(
        "--feature-types",
        default=",".join(DEFAULT_LIFTON_FEATURE_TYPES),
        help="Comma-separated parent feature types for LiftOn -f",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    feature_types = tuple(t.strip() for t in args.feature_types.split(",") if t.strip())

    lifton_config = LiftonRunConfig(
        ref_gff=args.ref_gff.resolve(),
        ref_fasta=args.ref_fasta.resolve(),
        target_fasta=args.target_fasta.resolve(),
        output_gff=args.output_gff.resolve(),
        threads=args.threads,
        executable=args.executable,
        feature_types=feature_types,
    )
    run_lifton(lifton_config)

    write_missing_gene_report(
        ref_gff=args.ref_gff.resolve(),
        projected_gff=args.output_gff.resolve(),
        output_report=args.missing_report.resolve(),
    )

    args.output_json.write_text(
        json.dumps(
            {
                "projected_gff": str(args.output_gff),
                "missing_report": str(args.missing_report),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
