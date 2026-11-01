#!/usr/bin/env python3
"""Sort BED12 models and copy the associated manifest."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-bed12", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    args = parser.parse_args()
    rows = [line for line in args.bed12.read_text().splitlines() if line.strip()]
    args.output_bed12.write_text("".join(f"{row}\n" for row in sorted(rows, key=lambda row: (row.split("\t")[0], int(row.split("\t")[1]), int(row.split("\t")[2]), row.split("\t")[3] if len(row.split("\t")) > 3 else ""))))
    shutil.copyfile(args.manifest, args.output_manifest)


if __name__ == "__main__":
    main()
