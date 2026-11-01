#!/usr/bin/env python3
"""Copy selected backend products into a deterministic publication directory."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--patterns", nargs="+", required=True)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    copied = set()
    for pattern in args.patterns:
        for product in sorted(args.source.glob(pattern)):
            if product.is_file() and product not in copied:
                shutil.copy2(product, args.destination / product.name)
                copied.add(product)


if __name__ == "__main__":
    main()
