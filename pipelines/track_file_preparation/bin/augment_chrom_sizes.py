#!/usr/bin/env python3
"""Append validated sequence lengths from a supplemental FASTA index."""

from __future__ import annotations

import argparse
import csv
import platform
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


def read_sizes(path: Path, fai: bool) -> List[Tuple[str, int]]:
    rows: List[Tuple[str, int]] = []
    seen = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            minimum_columns = 5 if fai else 2
            if len(fields) < minimum_columns:
                raise ValueError(f"Invalid sequence-size row at {path}:{line_number}: {line.rstrip()}")
            name, raw_length = fields[0], fields[1]
            if not name or not raw_length.isdigit() or int(raw_length) < 1:
                raise ValueError(f"Invalid sequence-size row at {path}:{line_number}: {line.rstrip()}")
            if name in seen:
                raise ValueError(f"Duplicate sequence name in {path}: {name}")
            seen.add(name)
            rows.append((name, int(raw_length)))
    if not rows:
        raise ValueError(f"No sequence-size rows found in {path}")
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--supplemental-fai", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--versions", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        primary = read_sizes(args.chrom_sizes, fai=False)
        supplemental = read_sizes(args.supplemental_fai, fai=True)
        combined: List[Tuple[str, int]] = list(primary)
        by_name: Dict[str, int] = dict(primary)
        report: List[List[str]] = [["sequence_name", "length", "source", "status"]]
        report.extend([[name, str(length), "chrom_sizes", "present"] for name, length in primary])
        for name, length in supplemental:
            if name in by_name:
                if by_name[name] != length:
                    raise ValueError(
                        f"Supplemental sequence {name} length {length} conflicts with primary chrom.sizes length {by_name[name]}"
                    )
                report.append([name, str(length), "supplemental_fasta", "already_present"])
                continue
            by_name[name] = length
            combined.append((name, length))
            report.append([name, str(length), "supplemental_fasta", "added"])

        with args.output.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(combined)
        with args.report.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(report)
        args.versions.write_text(
            '"AUGMENT_CHROM_SIZES":\n'
            f"    python: {platform.python_version()}\n"
        )
    except (OSError, ValueError, csv.Error) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
