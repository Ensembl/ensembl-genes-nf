#!/usr/bin/env python3
"""Write a deterministic identity record for a staged input file."""

import argparse
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    args.output.write_text(json.dumps({"label": args.label, "path": args.input.name, "sha256": digest}, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
