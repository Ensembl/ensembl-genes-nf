#!/usr/bin/env python3
"""Copy a GTF/status pair and write a backend audit."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def copy_if_needed(source: Path, destination: Path) -> None:
    if destination.exists() and source.samefile(destination):
        return
    shutil.copyfile(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--output-gtf", required=True, type=Path)
    parser.add_argument("--output-status", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()
    if not args.gtf.is_file() or args.gtf.stat().st_size == 0:
        raise SystemExit(f"{args.backend} GTF is missing or empty: {args.gtf}")
    copy_if_needed(args.gtf, args.output_gtf)
    copy_if_needed(args.status, args.output_status)
    status_rows = [line for line in args.status.read_text().splitlines() if line.strip()]
    state = status_rows[-1].split("\t")[-1] if status_rows else "EMPTY"
    args.audit.write_text(f"backend\tstatus\n{args.backend}\t{state}\n")


if __name__ == "__main__":
    main()
