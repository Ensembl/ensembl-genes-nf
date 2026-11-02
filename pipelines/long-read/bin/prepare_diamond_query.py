#!/usr/bin/env python3
"""Normalize a protein FASTA and classify whether it contains query records."""

import argparse
from pathlib import Path


def read_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    sequence: list[str] = []

    def flush() -> None:
        if header is not None and sequence:
            records.append((header, "".join(sequence)))

    for line in path.read_text().splitlines():
        value = line.strip()
        if not value:
            continue
        if value.startswith(">"):
            flush()
            header = value
            sequence = []
        elif header is not None:
            sequence.append(value)
    flush()
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--status", required=True, type=Path)
    args = parser.parse_args()

    records = read_records(args.input)
    with args.output.open("w") as output:
        for header, sequence in records:
            output.write(f"{header}\n{sequence}\n")
    state = "VALID" if records else "EMPTY"
    args.status.write_text(f"status\tstate\nquery\t{state}\n")


if __name__ == "__main__":
    main()
