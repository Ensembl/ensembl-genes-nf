#!/usr/bin/env python3
"""Extract splice junctions from a BAM, including clipped/indel-containing reads."""

import argparse
from collections import Counter


# pysam CIGAR operation codes.
REFERENCE_CONSUMING = {0, 2, 3, 7, 8}  # M, D, N, =, X
JUNCTION = 3               # N


def junctions_from_cigar(reference_start, cigartuples):
    """Return zero-based, half-open introns from a BAM CIGAR tuple list."""
    reference_position = reference_start
    junctions = []
    for operation, length in cigartuples or ():
        if operation == JUNCTION:
            junctions.append((reference_position, reference_position + length))
        if operation in REFERENCE_CONSUMING:
            reference_position += length
    return junctions


def extract_junctions(bam_path):
    """Return junction support keyed by contig, coordinates, and strand."""
    import pysam

    support = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue
            strand = "-" if read.is_reverse else "+"
            # BED chromosome names must match the BAM header exactly. FLAIR's
            # bundled junctions_from_sam helper adds "chr", but these inputs
            # use Ensembl-style contigs such as "1" and "MT".
            contig = read.reference_name
            for start, end in junctions_from_cigar(
                read.reference_start, read.cigartuples
            ):
                support[(contig, start, end, strand)] += 1
    return support


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("output")
    args = parser.parse_args()

    support = extract_junctions(args.bam)
    with open(args.output, "w") as output:
        for (contig, start, end, strand), count in sorted(support.items()):
            output.write(f"{contig}\t{start}\t{end}\t.\t{count}\t{strand}\n")


if __name__ == "__main__":
    main()
