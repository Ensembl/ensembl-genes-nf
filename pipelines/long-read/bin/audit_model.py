#!/usr/bin/env python3
"""Create a normalized model audit summary."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--validation-report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    models = "0"
    for line in args.validation_report.read_text().splitlines():
        fields = line.split("\t")
        if fields and fields[0] == "models" and len(fields) > 1:
            models = fields[1]
            break
    args.output.write_text(f"backend\tscope\tmodels\n{args.backend}\t{args.scope}\t{models}\n")


if __name__ == "__main__":
    main()
