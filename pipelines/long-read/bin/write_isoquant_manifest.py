#!/usr/bin/env python3
"""Write reproducible IsoQuant input/provenance metadata."""
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
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--genedb", type=Path)
    parser.add_argument("--data-type", required=True)
    parser.add_argument("--analysis", nargs="*", default=[])
    parser.add_argument("--large-output", nargs="*", default=[])
    parser.add_argument("--args", default="")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    inputs = {
        "reference": {"path": str(args.reference), "sha256": digest(args.reference)},
        "bam": {"path": str(args.bam), "sha256": digest(args.bam)},
        "data_type": args.data_type,
        "analysis": args.analysis,
        "large_output": args.large_output,
        "args": args.args,
    }
    if args.genedb:
        inputs["genedb"] = {"path": str(args.genedb), "sha256": digest(args.genedb)}
    args.output.write_text(json.dumps(inputs, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
