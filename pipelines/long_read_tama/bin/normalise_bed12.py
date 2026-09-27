#!/usr/bin/env python3
"""Normalise BED12 rows for downstream sequence extraction.

TAMA can emit zero-length terminal blocks.  They are not valid BED12 blocks
and bedtools getfasta rejects them, but removing them does not alter the
transcript span or any non-empty exon interval.
"""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output")
    parser.add_argument("report")
    args = parser.parse_args()

    rows = 0
    zero_length_blocks = 0
    with open(args.source) as source, open(args.output, "w") as output:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        for line_number, fields in enumerate(reader, 1):
            if not fields:
                continue
            if len(fields) != 12:
                raise SystemExit(f"Invalid BED12 field count on line {line_number}")
            block_count = int(fields[9])
            sizes = [int(value) for value in fields[10].rstrip(",").split(",")]
            starts = [int(value) for value in fields[11].rstrip(",").split(",")]
            if len(sizes) != block_count or len(starts) != block_count:
                raise SystemExit(f"BED12 block count mismatch on line {line_number}")

            kept = [(size, start) for size, start in zip(sizes, starts) if size > 0]
            zero_length_blocks += len(sizes) - len(kept)
            if not kept:
                raise SystemExit(f"BED12 row has no non-empty blocks on line {line_number}")

            span = int(fields[2]) - int(fields[1])
            for size, start in kept:
                if start < 0 or start + size > span:
                    raise SystemExit(f"BED12 block exceeds transcript span on line {line_number}")

            fields[9] = str(len(kept))
            fields[10] = ",".join(str(size) for size, _ in kept) + ","
            fields[11] = ",".join(str(start) for _, start in kept) + ","
            writer.writerow(fields)
            rows += 1

    with open(args.report, "w") as report:
        report.write("metric\tvalue\n")
        report.write(f"models\t{rows}\n")
        report.write(f"zero_length_blocks_removed\t{zero_length_blocks}\n")


if __name__ == "__main__":
    main()
