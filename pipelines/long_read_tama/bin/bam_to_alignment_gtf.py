#!/usr/bin/env python3
"""Convert primary mapped BAM alignments to exon-only, read-level GTF."""

import argparse
import re
import subprocess
from collections import defaultdict


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REFERENCE_CONSUMING = {"M", "D", "N", "=", "X"}
EXON_CONSUMING = {"M", "D", "=", "X"}


def alignment_exons(start, cigar):
    ref = start
    exon_start = ref
    exons = []
    for length_text, operation in CIGAR.findall(cigar):
        length = int(length_text)
        if operation == "N":
            if exon_start < ref:
                exons.append((exon_start, ref))
            ref += length
            exon_start = ref
        elif operation in EXON_CONSUMING:
            ref += length
        elif operation in REFERENCE_CONSUMING:
            ref += length
            if operation == "N":
                exon_start = ref
    if exon_start < ref:
        exons.append((exon_start, ref))
    return exons


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("output")
    args = parser.parse_args()

    seen = defaultdict(int)
    records = []
    command = ["samtools", "view", "-F", "2308", args.bam]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        read_id, flag_text, contig, pos_text, _mapq, cigar = fields[:6]
        if contig == "*" or cigar == "*":
            continue
        flag = int(flag_text)
        seen[read_id] += 1
        transcript_id = f"{read_id}.{seen[read_id]}"
        strand = "-" if flag & 16 else "+"
        for exon_number, (start, end) in enumerate(
            alignment_exons(int(pos_text), cigar), 1
        ):
            attrs = f'gene_id "{transcript_id}"; transcript_id "{transcript_id}"; exon_number "{exon_number}"; read_id "{read_id}";'
            records.append((contig, start, end, strand, attrs))

    records.sort(key=lambda row: (row[0], row[1], row[2], row[3], row[4]))
    with open(args.output, "w") as output:
        for contig, start, end, strand, attrs in records:
            output.write(f"{contig}\ttmerge\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")


if __name__ == "__main__":
    main()
