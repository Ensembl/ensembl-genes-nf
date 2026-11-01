#!/usr/bin/env python3
"""Write a JSON value supplied by a Nextflow value channel."""

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--value", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        value = json.loads(args.value)
    except json.JSONDecodeError as error:
        raise SystemExit(f"invalid JSON value: {error}") from error
    args.output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
