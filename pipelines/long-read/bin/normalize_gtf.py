#!/usr/bin/env python3
"""Normalize an optional GTF product and preserve its explicit status."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def has_exon(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    return any(
        len(fields := line.split("\t", 2)) == 3 and fields[2].split("\t", 1)[0] == "exon"
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--output-gtf", required=True, type=Path)
    parser.add_argument("--output-status", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--empty-message", required=True)
    args = parser.parse_args()
    if has_exon(args.gtf):
        shutil.copyfile(args.gtf, args.output_gtf)
        shutil.copyfile(args.status, args.output_status)
    else:
        args.output_gtf.write_text(f"# EMPTY: {args.empty_message}\n")
        args.output_status.write_text(f"{args.backend}\tEMPTY\n")


if __name__ == "__main__":
    main()
