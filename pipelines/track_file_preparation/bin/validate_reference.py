#!/usr/bin/env python3
"""Validate chromosome-size rows and write the reference check record."""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom-sizes", required=True, type=Path)
    parser.add_argument("--assembly-release", required=True)
    args = parser.parse_args()
    try:
        rows = args.chrom_sizes.read_text().splitlines()
        if not rows:
            raise ValueError("chrom.sizes is empty")
        for row in rows:
            fields = row.split("\t")
            if len(fields) != 2 or not fields[1].isdigit():
                raise ValueError(f"Invalid chrom.sizes row: {row}")
        Path("reference_check.tsv").write_text(
            "assembly_release\tchrom_sizes\tstatus\n"
            f"{args.assembly_release}\t{args.chrom_sizes}\tcomplete\n"
        )
        Path("reference.versions.yml").write_text(
            '"VALIDATE_REFERENCE":\n'
            f"    python: {platform.python_version()}\n"
        )
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
