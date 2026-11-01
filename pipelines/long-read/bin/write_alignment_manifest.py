#!/usr/bin/env python3
"""Write an alignment manifest with content identities for all inputs."""
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
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--bai", required=True, type=Path)
    parser.add_argument("--preset")
    parser.add_argument("--secondary", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    record = {
        "run_id": args.run_id,
        "input": args.input.name,
        "reference": args.reference.name,
        "reference_sha256": digest(args.reference),
        "input_sha256": digest(args.input),
        "bam_sha256": digest(args.bam),
        "bai_sha256": digest(args.bai),
    }
    if args.preset:
        record["minimap2_preset"] = args.preset
        record["secondary_mode"] = args.secondary
    args.output.write_text(json.dumps(record, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
