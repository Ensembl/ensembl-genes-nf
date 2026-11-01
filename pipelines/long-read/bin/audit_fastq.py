#!/usr/bin/env python3
"""Write a deterministic FASTQ audit record."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()
    if not args.reads.is_file() or args.reads.stat().st_size == 0:
        raise SystemExit(f"validated FASTQ is missing or empty: {args.reads}")
    args.audit.write_text(f"id\tstatus\tbytes\n{args.run_id}\tVALID\t{args.reads.stat().st_size}\n")


if __name__ == "__main__":
    main()
