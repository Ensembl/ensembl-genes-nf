#!/usr/bin/env python3
"""Extract deterministic read endpoints, splice junctions, and read features."""
from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
EXCLUDED_FLAGS = 4 + 256 + 512 + 1024 + 2048


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", type=Path)
    source.add_argument("--sam", type=Path)
    parser.add_argument("--read-features", required=True, type=Path)
    parser.add_argument("--five-prime", required=True, type=Path)
    parser.add_argument("--three-prime", required=True, type=Path)
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if args.bam is not None and shutil.which("samtools") is None:
        raise SystemExit("samtools is required for alignment feature extraction")
    if args.sam is not None:
        sam_text = args.sam.read_text()
    else:
        result = subprocess.run(
            ["samtools", "view", "-F", str(EXCLUDED_FLAGS), str(args.bam)],
            check=True,
            text=True,
            capture_output=True,
        )
        sam_text = result.stdout
    reads = []
    junctions: Counter[tuple[str, int, int, str]] = Counter()
    for line in sam_text.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        read_id, flag, chrom, position, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]
        strand = "-" if flag & 16 else "+"
        reference_position = position
        junction_count = 0
        for length_text, operation in CIGAR.findall(cigar):
            length = int(length_text)
            if operation == "N":
                junctions[(chrom, reference_position, reference_position + length, strand)] += 1
                junction_count += 1
                reference_position += length
            elif operation in "MD=X":
                reference_position += length
        end = reference_position
        aligned_bases = sum(int(length) for length, operation in CIGAR.findall(cigar) if operation in "M=X")
        five_prime = position if strand == "+" else end - 1
        three_prime = end - 1 if strand == "+" else position
        reads.append((read_id, chrom, strand, position, end, aligned_bases, junction_count, cigar, five_prime, three_prime))
    reads.sort(key=lambda row: (row[1], row[3], row[0]))
    with gzip.open(args.read_features, "wt") as handle:
        handle.write("read_id\tcontig\tstrand\tstart\tend\taligned_bases\tjunction_count\tcigar\n")
        for row in reads:
            handle.write("\t".join(map(str, row[:8])) + "\n")
    with args.five_prime.open("w") as five, args.three_prime.open("w") as three:
        for read_id, chrom, strand, _start, _end, _aligned, _junctions, _cigar, five_prime, three_prime in reads:
            five.write(f"{chrom}\t{five_prime}\t{five_prime + 1}\t{read_id}\t0\t{strand}\n")
            three.write(f"{chrom}\t{three_prime}\t{three_prime + 1}\t{read_id}\t0\t{strand}\n")
    with args.junctions.open("w") as handle:
        for (chrom, start, end, strand), support in sorted(junctions.items()):
            name = f"{chrom}:{start}-{end}"
            handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{support}\t{strand}\n")
    args.summary.write_text(
        "run_id\treads\tjunctions\n"
        f"{args.run_id}\t{len(reads)}\t{len(junctions)}\n"
    )


if __name__ == "__main__":
    main()
