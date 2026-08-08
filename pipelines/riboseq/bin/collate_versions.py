#!/usr/bin/env python3
"""Merge per-task Nextflow versions YAML files into one deterministic report."""

from __future__ import annotations

import sys
from pathlib import Path

import yaml


def main() -> int:
    merged: dict[str, dict[str, str]] = {}
    for raw_path in sys.argv[1:]:
        path = Path(raw_path)
        for data in yaml.safe_load_all(path.read_text()):
            data = data or {}
            for process, tools in data.items():
                merged.setdefault(str(process), {}).update(
                    {str(tool): str(version) for tool, version in (tools or {}).items()}
                )

    for process in sorted(merged):
        print(f"{process}:")
        for tool in sorted(merged[process]):
            print(f"  {tool}: {merged[process][tool]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
