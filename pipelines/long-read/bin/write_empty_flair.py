#!/usr/bin/env python3
"""Write the explicit empty FLAIR result for a junction-free shard."""

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--junction-status", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--shard-id", required=True)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--products", required=True, type=Path)
    args = parser.parse_args()
    rows = [line for line in args.junction_status.read_text().splitlines() if line.strip()]
    if not rows or rows[-1].split("\t")[-1] != "EMPTY":
        raise SystemExit(f"expected EMPTY junction status: {args.junction_status}")
    args.gtf.write_text(f"# EMPTY: no splice junctions for shard {args.shard_id}\n")
    args.status.write_text(f"flair\tEMPTY\t{args.shard_id}\n")
    args.products.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps({
        "backend": "flair",
        "run_id": args.run_id,
        "shard_id": args.shard_id,
        "status": "EMPTY",
        "inputs": [],
        "native_products": [],
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
