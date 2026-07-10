#!/usr/bin/env python3
"""Create a ChOROS BAM with conservative transcript assignments and ZW weights.

Input must be query-name grouped (for example, ``samtools collate`` output).
Only fragments with one primary transcript alignment are retained.  Collapsed
read abundance encoded as a trailing ``_xN`` is written to the ZW:f tag that
ChOROS understands.
"""

from __future__ import annotations

import argparse
import re
from itertools import groupby

import pysam


COUNT_RE = re.compile(r"_x(\d+)$")


def collapsed_count(query_name: str) -> int:
    match = COUNT_RE.search(query_name)
    return int(match.group(1)) if match else 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metrics", required=True)
    args = parser.parse_args()

    total = retained = ambiguous = 0
    with pysam.AlignmentFile(args.bam, "rb") as source:
        with pysam.AlignmentFile(args.output, "wb", template=source) as dest:
            grouped = groupby(source.fetch(until_eof=True), key=lambda read: read.query_name)
            for _, records_iter in grouped:
                total += 1
                records = [
                    read
                    for read in records_iter
                    if not read.is_unmapped and not read.is_secondary and not read.is_supplementary
                ]
                if len(records) != 1:
                    ambiguous += 1
                    continue
                read = records[0]
                read.set_tag("ZW", float(collapsed_count(read.query_name)), value_type="f")
                dest.write(read)
                retained += 1

    with open(args.metrics, "w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"query_groups_total\t{total}\n")
        handle.write(f"query_groups_retained\t{retained}\n")
        handle.write(f"query_groups_ambiguous\t{ambiguous}\n")
        fraction = retained / total if total else 0.0
        handle.write(f"retained_fraction\t{fraction:.8f}\n")


if __name__ == "__main__":
    main()
