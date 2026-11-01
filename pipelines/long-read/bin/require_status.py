#!/usr/bin/env python3
"""Require a status record containing a requested value."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--value", required=True)
    args = parser.parse_args()
    if args.value not in args.status.read_text().split():
        raise SystemExit(f"{args.status} does not contain required status {args.value}")


if __name__ == "__main__":
    main()
