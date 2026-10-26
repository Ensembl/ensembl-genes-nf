#!/usr/bin/env python3
"""Copy version-dependent FLAIR products to stable comparison names."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def find_one(
    root: Path,
    suffixes: tuple[str, ...],
    prefix: str | None = None,
    exclude: Path | None = None,
) -> Path:
    excluded = exclude.resolve() if exclude is not None else None
    matches = sorted(
        {
            path
            for suffix in suffixes
            for path in root.rglob(f"*{suffix}")
            if path.is_file()
            and (prefix is None or path.name.startswith(prefix))
            and (excluded is None or excluded not in path.resolve().parents)
        }
    )
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one FLAIR product with suffixes {suffixes}, found {matches}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("prefix", nargs="?", default=None)
    args = parser.parse_args()
    # FLAIR may create this directory before the stable-name copy runs. Exclude
    # it from discovery so prior normalized outputs cannot be mistaken for inputs.
    products = {
        "isoforms.bed": find_one(args.source, (".isoforms.bed", ".bed"), args.prefix, args.destination),
        "isoforms.gtf": find_one(args.source, (".isoforms.gtf", ".gtf"), args.prefix, args.destination),
        "isoforms.fa": find_one(args.source, (".isoforms.fa", ".fa"), args.prefix, args.destination),
        "read.map.txt": find_one(
            args.source, (".read.map.txt", ".isoform.map.txt"), args.prefix, args.destination
        ),
    }
    args.destination.mkdir(parents=True, exist_ok=True)
    for name, source in products.items():
        shutil.copy2(source, args.destination / name)
    (args.destination / "product_manifest.tsv").write_text(
        "stable_name\tsource_name\n"
        + "".join(f"{name}\t{source.relative_to(args.source)}\n" for name, source in products.items())
    )


if __name__ == "__main__":
    main()
