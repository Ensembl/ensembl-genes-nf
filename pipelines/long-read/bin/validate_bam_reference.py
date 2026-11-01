#!/usr/bin/env python3
"""Verify that BAM contigs and lengths exactly match the reference dictionary."""
from __future__ import annotations

import argparse
from pathlib import Path


def read_dictionary(path: Path) -> list[tuple[str, int]]:
    return [(fields[0], int(fields[1])) for fields in (line.split("\t") for line in path.read_text().splitlines() if line.strip())]


def read_idxstats(path: Path) -> list[tuple[str, int]]:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[0] != "*":
            rows.append((fields[0], int(fields[1])))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--idxstats", required=True, type=Path)
    parser.add_argument("--reference-dictionary", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    reference = read_dictionary(args.reference_dictionary)
    bam = read_idxstats(args.idxstats)
    errors = []
    if reference != bam:
        errors.append(f"reference/BAM dictionaries differ: reference={reference!r} bam={bam!r}")
    args.report.write_text("metric\tvalue\n" + f"contigs\t{len(bam)}\nerrors\t{len(errors)}\n" + "".join(f"error\t{error}\n" for error in errors))
    if errors:
        raise SystemExit(errors[0])


if __name__ == "__main__":
    main()
