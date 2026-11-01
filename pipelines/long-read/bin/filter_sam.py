#!/usr/bin/env python3
"""Remove program-record headers from a SAM stream while preserving records."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.input.open() as source, args.output.open("w") as target:
        for line in source:
            if not line.startswith("@PG\t"):
                target.write(line)


if __name__ == "__main__":
    main()
