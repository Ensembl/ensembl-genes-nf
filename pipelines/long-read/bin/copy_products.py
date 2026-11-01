#!/usr/bin/env python3
"""Copy a declared list of optional products into a destination directory."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="*", type=Path, default=[])
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    for product in args.inputs:
        if product.is_file():
            shutil.copy2(product, args.destination / product.name)


if __name__ == "__main__":
    main()
