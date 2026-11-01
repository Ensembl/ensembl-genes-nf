#!/usr/bin/env python3
"""Create placeholder transcript FASTA records from BED12 model IDs."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.bed12.open() as source, args.output.open("w") as output:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3]:
                output.write(f">{fields[3]}\nNNNN\n")


if __name__ == "__main__":
    main()
