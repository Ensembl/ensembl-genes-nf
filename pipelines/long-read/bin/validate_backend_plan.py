#!/usr/bin/env python3
"""Validate an audit report and write a structured result."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not args.audit.is_file() or args.audit.stat().st_size == 0:
        raise SystemExit(f"backend audit is missing or empty: {args.audit}")
    errors = sum(1 for line in args.audit.read_text().splitlines() if line.startswith("# error\t"))
    args.output.write_text(f"status\tvalue\nVALID\t{'true' if errors == 0 else 'false'}\nerrors\t{errors}\n")
    if errors:
        raise SystemExit(f"backend audit contains {errors} errors")


if __name__ == "__main__":
    main()
