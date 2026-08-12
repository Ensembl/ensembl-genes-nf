#!/usr/bin/env python3
"""Convert legacy pairwise metric headers to the public source A/source B contract."""

from __future__ import annotations

import argparse
from pathlib import Path


REPLACEMENTS = (
    ("ens_to_cat", "source_a_to_source_b"),
    ("cat_to_ens", "source_b_to_source_a"),
    ("ensembl", "source_a"),
    ("cat", "source_b"),
)


def normalise_header(header: str) -> str:
    result = header
    for old, new in REPLACEMENTS:
        result = result.replace(old, new)
    return result


def normalise(input_path: Path, output_path: Path) -> None:
    with input_path.open(encoding="utf-8") as source, output_path.open(
        "w", encoding="utf-8"
    ) as target:
        first = True
        for line in source:
            if first and line.strip():
                target.write(normalise_header(line))
                first = False
            else:
                target.write(line)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    normalise(args.input, args.output)


if __name__ == "__main__":
    main()
