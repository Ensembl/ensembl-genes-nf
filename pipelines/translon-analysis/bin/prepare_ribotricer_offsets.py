#!/usr/bin/env python3
"""Convert a RiboMetric/RiboWaltz offset table to Ribotricer's CLI lists."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_offsets(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs: list[tuple[int, int]] = []
    for row in rows:
        length = row.get("length") or row.get("read_length") or row.get("read_len")
        offset = row.get("offset") or row.get("new_offset") or row.get("computed_offset")
        if not length or not offset:
            continue
        try:
            pairs.append((int(float(length)), int(float(offset))))
        except ValueError:
            continue
    pairs = sorted(set(pairs))
    if not pairs:
        raise ValueError(f"No read-length/offset pairs found in {path}")
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--read-lengths", required=True, type=Path)
    parser.add_argument("--offsets", required=True, type=Path)
    args = parser.parse_args()
    pairs = read_offsets(args.input)
    args.read_lengths.write_text(",".join(str(length) for length, _ in pairs) + "\n")
    args.offsets.write_text(",".join(str(offset) for _, offset in pairs) + "\n")


if __name__ == "__main__":
    main()
