#!/usr/bin/env python3
"""Select a backend GTF product, write its explicit status, and handle empties."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def has_exon(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    return any(
        len(fields := line.rstrip("\n").split("\t")) >= 3 and fields[2] == "exon"
        for line in path.open()
        if line and not line.startswith("#")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--source", type=Path)
    source.add_argument("--source-dir", type=Path)
    parser.add_argument("--pattern")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--empty-message", required=True)
    parser.add_argument("--require-source", action="store_true")
    args = parser.parse_args()
    source_path = args.source
    if args.source_dir:
        if not args.pattern:
            raise SystemExit("--pattern is required with --source-dir")
        source_path = next(iter(sorted(args.source_dir.rglob(args.pattern))), None)
    if source_path is None:
        raise SystemExit("no GTF source was selected")
    if args.require_source and not source_path.is_file():
        raise SystemExit(f"required GTF product is missing: {source_path}")
    if source_path.is_file() and has_exon(source_path):
        shutil.copyfile(source_path, args.output)
        state = "SUCCESS"
    else:
        args.output.write_text(f"# EMPTY: {args.empty_message}\n")
        state = "EMPTY"
    args.status.write_text(f"{args.backend}\t{state}\n")


if __name__ == "__main__":
    main()
