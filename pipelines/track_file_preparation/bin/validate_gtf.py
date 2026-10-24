#!/usr/bin/env python3
"""Validate the structural fields used by the track-preparation GTF path."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, List


def validate_gtf(path: Path) -> List[str]:
    errors: List[str] = []
    data_rows = 0

    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\r\n")
            if not stripped.strip() or stripped.lstrip().startswith("#"):
                continue

            data_rows += 1
            fields = stripped.split("\t")
            if len(fields) != 9:
                errors.append(
                    f"line {line_number}: expected 9 tab-separated columns, found {len(fields)}"
                )
                continue

            if not all(fields[index] for index in (0, 1, 2, 5, 6, 7, 8)):
                errors.append(f"line {line_number}: required GTF field is empty")
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                errors.append(f"line {line_number}: start and end must be integers")
                continue

            if start < 1 or end < start:
                errors.append(
                    f"line {line_number}: coordinates must satisfy 1 <= start <= end"
                )

    if data_rows == 0:
        errors.append("GTF has no data rows")
    return errors


def main(argv: Iterable[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 1:
        print(f"Usage: {Path(sys.argv[0]).name} GTF", file=sys.stderr)
        return 2

    path = Path(args[0])
    try:
        errors = validate_gtf(path)
    except OSError as exc:
        print(f"Cannot read GTF {path}: {exc}", file=sys.stderr)
        return 2

    for error in errors:
        print(f"{path}: {error}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
