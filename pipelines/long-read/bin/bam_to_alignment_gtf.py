#!/usr/bin/env python3
"""Convert a name-sorted/read-alignment BAM through samtools view records.

The production container supplies pysam or samtools; this helper deliberately
fails when neither is available instead of manufacturing a biological model.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import re
from pathlib import Path

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
EXCLUDED_FLAGS = 4 + 256 + 512 + 1024 + 2048


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", type=Path)
    source.add_argument("--sam", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.bam is not None:
        if shutil.which("samtools") is None:
            raise SystemExit("--bam requires samtools; use --sam in the Python-only process")
        source_handle = subprocess.Popen(["samtools", "view", "-F", str(EXCLUDED_FLAGS), str(args.bam)], stdout=subprocess.PIPE, text=True).stdout
    else:
        source_handle = args.sam.open()
    with args.output.open("w") as handle:
        for line in source_handle:
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            read_id, flag, chrom, start, cigar = fields[0], int(fields[1]), fields[2], int(fields[3]) - 1, fields[5]
            strand = "-" if flag & 16 else "+"
            reference_position = start
            exon_start = start
            exon_number = 0
            for length_text, operation in CIGAR.findall(cigar):
                length = int(length_text)
                if operation in "M=X":
                    reference_position += length
                elif operation == "N":
                    if reference_position > exon_start:
                        exon_number += 1
                        handle.write(f'{chrom}\tpipeline\texon\t{exon_start + 1}\t{reference_position}\t.\t{strand}\t.\tgene_id "{read_id}"; transcript_id "{read_id}"; exon_number "{exon_number}";\n')
                    reference_position += length
                    exon_start = reference_position
                elif operation in "D":
                    reference_position += length
            if reference_position > exon_start:
                exon_number += 1
                handle.write(f'{chrom}\tpipeline\texon\t{exon_start + 1}\t{reference_position}\t.\t{strand}\t.\tgene_id "{read_id}"; transcript_id "{read_id}"; exon_number "{exon_number}";\n')
            if exon_number == 0:
                continue
    if args.bam is not None:
        source_handle.close()
if __name__ == "__main__":
    main()
