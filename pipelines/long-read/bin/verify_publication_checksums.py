#!/usr/bin/env python3
"""Verify a relative-path SHA-256 publication checksum table."""

from __future__ import annotations

import argparse
import hashlib
import csv
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--checksums", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    errors = 0
    with args.checksums.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            # publishDir may copy asynchronously after this finalization task.
            # Verify the staged source while reporting the eventual published path.
            relative = Path(row.get("source_path") or row["path"])
            reported_path = row["path"]
            expected = row["sha256"]
            target = args.root / relative
            if not target.is_file():
                actual = "MISSING"
                status = "ERROR"
            else:
                actual = digest(target)
                status = "VALID" if actual == expected else "ERROR"
            errors += status == "ERROR"
            rows.append({"path": reported_path, "source_path": str(relative), "expected_sha256": expected, "actual_sha256": actual, "status": status})

    args.report.write_text(
        "path\tsource_path\texpected_sha256\tactual_sha256\tstatus\n"
        + "".join(f"{row['path']}\t{row['source_path']}\t{row['expected_sha256']}\t{row['actual_sha256']}\t{row['status']}\n" for row in rows)
        + f"# errors\t{errors}\n"
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
