#!/usr/bin/env python3
"""Copy version-dependent FLAIR products to stable comparison names."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def find_one(root: Path, suffixes: tuple[str, ...], prefix: str | None = None) -> Path:
    matches = sorted(
        {
            path
            for suffix in suffixes
            for path in root.rglob(f"*{suffix}")
            if path.is_file() and (prefix is None or path.name.startswith(prefix))
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
    args.destination.mkdir(parents=True, exist_ok=False)
    products = {
        "isoforms.bed": find_one(args.source, (".isoforms.bed", ".bed"), args.prefix),
        "isoforms.gtf": find_one(args.source, (".isoforms.gtf", ".gtf"), args.prefix),
        "isoforms.fa": find_one(args.source, (".isoforms.fa", ".fa"), args.prefix),
        "read.map.txt": find_one(args.source, (".read.map.txt", ".isoform.map.txt"), args.prefix),
    }
    for name, source in products.items():
        shutil.copy2(source, args.destination / name)
    (args.destination / "product_manifest.tsv").write_text(
        "stable_name\tsource_name\n"
        + "".join(f"{name}\t{source.relative_to(args.source)}\n" for name, source in products.items())
    )


if __name__ == "__main__":
    main()
