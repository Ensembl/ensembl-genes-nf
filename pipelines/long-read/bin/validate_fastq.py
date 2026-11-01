#!/usr/bin/env python3
"""Validate plain or gzip-compressed FASTQ without shell text parsing."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def records(path: Path):
    opener = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
    with opener(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError("truncated FASTQ record")
            yield header.rstrip("\r\n"), sequence.rstrip("\r\n"), plus.rstrip("\r\n"), quality.rstrip("\r\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--status", required=True, type=Path)
    args = parser.parse_args()
    if not args.reads.is_file() or args.reads.stat().st_size == 0:
        raise SystemExit(f"FASTQ is missing or empty: {args.reads}")
    count = 0
    try:
        for header, sequence, plus, quality in records(args.reads):
            count += 1
            if not header.startswith("@"):
                raise ValueError(f"record {count} does not start with @")
            if not plus.startswith("+"):
                raise ValueError(f"record {count} does not contain a + separator")
            if len(sequence) != len(quality):
                raise ValueError(f"record {count} sequence and quality lengths differ")
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile, ValueError) as error:
        raise SystemExit(f"invalid FASTQ {args.reads}: {error}") from error
    if count == 0:
        raise SystemExit(f"FASTQ contains no records: {args.reads}")
    args.status.write_text(f"id\tstatus\trecords\n{args.run_id}\tVALID\t{count}\n")


if __name__ == "__main__":
    main()
