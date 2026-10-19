#!/usr/bin/env python3
"""Extract primary-read intron chains and junction support from a sorted BAM."""
import sys
from collections import Counter
import pysam

def main():
    bam_path, sample, summary_path, chains_path, junctions_path = sys.argv[1:]
    chains = Counter()
    junctions = Counter()
    mapped = spliced = 0
    intron_total = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam, open(chains_path, "w") as chains_out:
        chains_out.write("read_id\tcontig\tstrand\tintron_count\tintron_chain\tmapq\n")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.reference_name:
                continue
            mapped += 1
            refpos = read.reference_start
            introns = []
            for operation, length in read.cigartuples or []:
                if operation == 3:
                    introns.append((refpos, refpos + length))
                    junctions[(read.reference_name, refpos, refpos + length, "-" if read.is_reverse else "+")] += 1
                if operation in {0, 2, 3, 7, 8}:
                    refpos += length
            if introns:
                spliced += 1
                intron_total += len(introns)
            chain = ";".join(f"{start}-{end}" for start, end in introns)
            chains[(read.reference_name, "-" if read.is_reverse else "+", chain)] += 1
            chains_out.write(f"{read.query_name}\t{read.reference_name}\t{'-' if read.is_reverse else '+'}\t{len(introns)}\t{chain}\t{read.mapping_quality}\n")
    with open(junctions_path, "w") as out:
        for (contig, start, end, strand), support in sorted(junctions.items()):
            out.write(f"{contig}\t{start}\t{end}\t.\t{support}\t{strand}\n")
    unique_chains = sum(1 for (_, _, chain) in chains if chain)
    with open(summary_path, "w") as out:
        out.write("sample_id\tprimary_mapped_reads\tspliced_reads\tunique_intron_chain_count\tunique_junction_count\tmean_introns_per_spliced_read\n")
        mean = intron_total / spliced if spliced else 0.0
        out.write(f"{sample}\t{mapped}\t{spliced}\t{unique_chains}\t{len(junctions)}\t{mean:.6f}\n")

if __name__ == "__main__":
    main()
