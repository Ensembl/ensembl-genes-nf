#!/usr/bin/env python3
"""Write a checksum-bearing Bambu input manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--track-reads", required=True, type=str)
    parser.add_argument("--quantify", required=True, type=str)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = {
        "backend": "bambu",
        "bam": {"path": str(args.bam), "sha256": digest(args.bam)},
        "reference": {"path": str(args.reference), "sha256": digest(args.reference)},
        "track_reads": args.track_reads.lower() == "true",
        "quantify": args.quantify.lower() == "true",
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
