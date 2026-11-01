#!/usr/bin/env python3
"""Extract primary passing-read splice junctions from a BAM via samtools."""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
EXCLUDED_FLAGS = 4 + 256 + 512 + 1024 + 2048


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if shutil.which("samtools") is None:
        raise SystemExit("samtools is required for BAM junction extraction")
    result = subprocess.run(
        ["samtools", "view", "-F", str(EXCLUDED_FLAGS), str(args.bam)],
        check=True,
        text=True,
        capture_output=True,
    )
    junctions: Counter[tuple[str, int, int, str]] = Counter()
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        chrom, position, cigar = fields[2], int(fields[3]) - 1, fields[5]
        strand = "-" if int(fields[1]) & 16 else "+"
        reference_position = position
        for length_text, operation in CIGAR.findall(cigar):
            length = int(length_text)
            if operation == "N":
                junctions[(chrom, reference_position, reference_position + length, strand)] += 1
                reference_position += length
            elif operation in "MD=X":
                reference_position += length

    with args.output.open("w") as handle:
        for (chrom, start, end, strand), support in sorted(junctions.items()):
            name = f"{chrom}:{start}-{end}"
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{support}\t{strand}\n")


if __name__ == "__main__":
    main()
