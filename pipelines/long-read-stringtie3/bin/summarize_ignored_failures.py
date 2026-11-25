#!/usr/bin/env python3
"""Summarize terminally failed (and ignored) tasks from a Nextflow trace."""

import argparse
import csv
from pathlib import Path


def attempt_number(row: dict[str, str]) -> int:
    try:
        return int(row.get("attempt", "1") or "1")
    except ValueError:
        return 1


def error_excerpt(workdir: str) -> str:
    if not workdir:
        return ""
    path = Path(workdir) / ".command.err"
    if not path.is_file():
        return ""
    return path.read_text(errors="replace")[-3000:].strip().replace("\n", " | ")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path, help="Nextflow trace file in TSV format")
    parser.add_argument("output", type=Path, help="Output TSV for terminal failures")
    args = parser.parse_args()

    with args.trace.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    latest: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row.get("process", ""), row.get("tag", ""))
        if key not in latest or attempt_number(row) >= attempt_number(latest[key]):
            latest[key] = row

    failures = [
        row for row in latest.values()
        if row.get("status", "").upper() == "FAILED"
    ]
    failures.sort(key=lambda row: (row.get("process", ""), row.get("tag", "")))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["process", "tag", "attempt", "status", "exit", "workdir", "error_excerpt"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in failures:
            writer.writerow({**row, "error_excerpt": error_excerpt(row.get("workdir", ""))})

    print(f"Wrote {len(failures)} terminal failed tasks to {args.output}")


if __name__ == "__main__":
    main()
