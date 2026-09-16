#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("bam")
    p.add_argument("bai")
    p.add_argument("output")
    p.add_argument("threshold", type=int)
    a = p.parse_args()
    if not Path(a.bam).is_file() or not Path(a.bai).is_file():
        raise SystemExit("BAM and BAI must be readable")
    subprocess.run(["samtools", "quickcheck", "-v", a.bam], check=True)
    result = subprocess.run(["samtools", "idxstats", a.bam], check=True, text=True, capture_output=True)
    lines = ["contig\treference_bases\tmapped_reads\teligible\tresource_class"]
    for raw in result.stdout.splitlines():
        fields = raw.split("\t")
        if len(fields) < 4:
            continue
        contig, length, mapped = fields[0], int(fields[1]), int(fields[2])
        if contig == "*" or length <= 0:
            continue
        resource = "very_large" if mapped >= a.threshold * 5 else "large" if mapped >= a.threshold else "small"
        lines.append(f"{contig}\t{length}\t{mapped}\ttrue\t{resource}")
    if len(lines) == 1:
        raise SystemExit("BAM has no eligible reference contigs")
    Path(a.output).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
