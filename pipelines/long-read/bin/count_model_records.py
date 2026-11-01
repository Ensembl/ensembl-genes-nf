#!/usr/bin/env python3
"""Write small line-count reports for model and evidence products."""
from __future__ import annotations

import argparse
from pathlib import Path


def count_records(path: Path, skip_header: bool = False) -> int:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if skip_header and lines and lines[0].lower().startswith(("query\t", "status\t")):
        return len(lines) - 1
    return len(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", type=Path)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--classification-policy", default="longest_atg_orf")
    parser.add_argument("--hits-only", action="store_true")
    args = parser.parse_args()
    if args.hits_only:
        args.report.write_text(f"status\tmodels\nCOMPLETE\t{count_records(args.hits, skip_header=True)}\n")
        return
    if args.models is None:
        raise SystemExit("--models is required unless --hits-only is used")
    models = count_records(args.models)
    hits = count_records(args.hits, skip_header=True)
    args.report.write_text(
        "backend\tscope\tmodels\tdiamond_hits\tclassification_policy\n"
        f"{args.backend}\t{args.scope}\t{models}\t{hits}\t{args.classification_policy}\n"
    )


if __name__ == "__main__":
    main()
