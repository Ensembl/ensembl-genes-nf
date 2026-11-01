#!/usr/bin/env python3
"""Create content-addressed IDs for evidence and assessment datasets."""
from __future__ import annotations

import argparse
import hashlib
import json


def dataset_id(kind: str, payload: dict) -> str:
    # Compatibility and matcher policy fields are intentionally part of the
    # digest: a dictionary-only assessment is not the same dataset as one that
    # also permits sequence matching.
    canonical = json.dumps({"kind": kind, **payload}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", required=True)
    parser.add_argument("--payload", required=True, help="JSON object")
    args = parser.parse_args()
    print(dataset_id(args.kind, json.loads(args.payload)))


if __name__ == "__main__":
    main()
