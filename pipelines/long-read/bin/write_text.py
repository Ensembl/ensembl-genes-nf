#!/usr/bin/env python3
"""Decode a base64 text value supplied by a Nextflow value channel."""

import argparse
import base64
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base64", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        value = base64.b64decode(args.base64, validate=True).decode()
    except (ValueError, UnicodeDecodeError) as error:
        raise SystemExit(f"invalid base64 text: {error}") from error
    args.output.write_text(value)


if __name__ == "__main__":
    main()
