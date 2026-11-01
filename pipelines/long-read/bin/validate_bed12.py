#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--reference-dictionary", type=Path)
    args = parser.parse_args()
    contigs = {}
    if args.reference_dictionary:
        for line in args.reference_dictionary.read_text().splitlines():
            name, length = line.split("\t")[:2]
            contigs[name] = int(length)
    errors = []
    count = 0
    with args.bed12.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                errors.append(f"line {number}: expected 12 columns")
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            block_count = int(fields[9])
            sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
            starts = [int(value) for value in fields[11].rstrip(",").split(",") if value]
            if start < 0 or end <= start or block_count != len(sizes) or block_count != len(starts):
                errors.append(f"line {number}: invalid coordinates or block lists")
                continue
            intervals = [(start + offset, start + offset + size) for offset, size in zip(starts, sizes)]
            if any(left < 0 or right <= left or right > end for left, right in intervals):
                errors.append(f"line {number}: block outside transcript bounds")
            if any(right >= next_left for (_left, right), (next_left, _next_right) in zip(intervals, intervals[1:])):
                errors.append(f"line {number}: overlapping or touching blocks")
            if fields[5] not in "+-.":
                errors.append(f"line {number}: invalid strand")
            if contigs and (chrom not in contigs or end > contigs[chrom]):
                errors.append(f"line {number}: contig missing or coordinate exceeds reference")
            count += 1
    with args.report.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"models\t{count}\n")
        handle.write(f"errors\t{len(errors)}\n")
        for error in errors:
            handle.write(f"error\t{error}\n")
    if errors:
        raise SystemExit("invalid BED12: " + "; ".join(errors[:3]))


if __name__ == "__main__":
    main()
