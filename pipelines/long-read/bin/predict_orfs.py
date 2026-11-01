#!/usr/bin/env python3
"""Apply the pipeline's deterministic placeholder ORF transformation."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.transcripts.open() as source, args.output.open("w") as output:
        for line in source:
            if line.startswith(">"):
                output.write(line)
            else:
                sequence = line.rstrip("\r\n")
                output.write(("M" + sequence[1:] if sequence else "") + "\n")


if __name__ == "__main__":
    main()
