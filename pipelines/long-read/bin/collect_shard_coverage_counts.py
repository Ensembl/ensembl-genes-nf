#!/usr/bin/env python3
"""Combine per-source and per-shard mapped-read counts."""

import argparse
from pathlib import Path


def read_count(path: Path) -> int:
    return int(path.read_text().strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-count", required=True, type=Path)
    parser.add_argument("--shard-counts", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source = read_count(args.source_count)
    shards = sum(read_count(path) for path in args.shard_counts)
    args.output.write_text(f"source_mapped_reads\t{source}\nshard_mapped_reads\t{shards}\n")


if __name__ == "__main__":
    main()
