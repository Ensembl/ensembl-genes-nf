#!/usr/bin/env python3
"""Concatenate GTF products while preserving an explicit empty state."""
from __future__ import annotations

import argparse
from pathlib import Path


def has_exon(path: Path) -> bool:
    return any(
        len(fields := line.rstrip("\n").split("\t")) >= 3 and fields[2] == "exon"
        for line in path.read_text().splitlines()
        if line and not line.startswith("#")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--scope")
    parser.add_argument("--skipped", type=Path)
    parser.add_argument("--merge-status", type=Path)
    parser.add_argument("--sort", choices=["name", "coordinates"], default="name")
    args = parser.parse_args()
    records = []
    paths = sorted(args.inputs, key=lambda item: item.name)
    for path in paths:
        if path.is_file():
            records.extend(line for line in path.read_text().splitlines() if line and not line.startswith("#"))
    if args.sort == "coordinates":
        records.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[3]), int(line.split("\t")[4])))
    if records:
        args.output.write_text("".join(f"{line}\n" for line in records))
        state = "COMPLETE"
    else:
        args.output.write_text(f"# EMPTY: no native exon records for {args.backend}\n")
        state = "EMPTY"
    args.status.write_text(f"backend\tstatus\n{args.backend}\t{state}\n")
    if args.skipped:
        skipped = "COMPLETE" if state == "COMPLETE" else "SKIPPED_EMPTY"
        args.skipped.write_text(f"input\tstatus\n{args.backend}\t{skipped}\n")
    if args.merge_status:
        if not args.scope:
            raise SystemExit("--merge-status requires --scope")
        args.merge_status.write_text(f"backend\tscope\tstatus\n{args.backend}\t{args.scope}\t{state}\n")


if __name__ == "__main__":
    main()
