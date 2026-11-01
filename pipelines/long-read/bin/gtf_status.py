#!/usr/bin/env python3
"""Classify a GTF-like product by the presence of exon records."""
from __future__ import annotations

import argparse
from pathlib import Path


def has_exon(path: Path) -> bool:
    if not path.is_file():
        return False
    for line in path.read_text().splitlines():
        if line and not line.startswith("#") and len(line.split("\t")) >= 3 and line.split("\t", 2)[2] == "exon":
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--empty-message", required=True)
    args = parser.parse_args()
    state = "SUCCESS" if has_exon(args.gtf) else "EMPTY"
    if state == "EMPTY":
        args.gtf.write_text(f"# EMPTY: {args.empty_message}\n")
    args.status.write_text(f"backend\tstatus\n{args.backend}\t{state}\n")


if __name__ == "__main__":
    main()
