#!/usr/bin/env python3
"""Summarise primary-alignment identity from minimap2 long cs tags."""
import re
import statistics
import sys
import pysam

CS_TOKEN = re.compile(r"(:\d+|\*[a-z][a-z]|\+[A-Za-z]+|-[A-Za-z]+)")

def counts(cs):
    values = [0, 0, 0, 0]
    for token in CS_TOKEN.findall(cs or ""):
        if token.startswith(":"):
            values[0] += int(token[1:])
        elif token.startswith("*"):
            values[1] += 1
        elif token.startswith("+"):
            values[2] += len(token) - 1
        elif token.startswith("-"):
            values[3] += len(token) - 1
    return values

def main():
    bam_path, sample, summary_path, reads_path = sys.argv[1:]
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            matches, substitutions, inserted, deleted = counts(read.get_tag("cs")) if read.has_tag("cs") else (0, 0, 0, 0)
            denominator = matches + substitutions + inserted + deleted
            identity = matches / denominator if denominator else 0.0
            aligned_query_bases = matches + substitutions + inserted
            rows.append((identity, matches, substitutions, inserted, deleted, aligned_query_bases, read.mapping_quality, read.query_name, read.reference_name or "*"))
    with open(reads_path, "w") as out:
        out.write("read_id\tcontig\tidentity\tmatches\tsubstitutions\tinserted_bases\tdeleted_bases\taligned_query_bases\tmapq\n")
        for identity, matches, substitutions, inserted, deleted, length, mapq, name, contig in rows:
            out.write(f"{name}\t{contig}\t{identity:.8f}\t{matches}\t{substitutions}\t{inserted}\t{deleted}\t{length}\t{mapq}\n")
    identities = [row[0] for row in rows]
    weighted_denominator = sum(row[1] + row[2] + row[3] + row[4] for row in rows)
    weighted_mean = sum(row[1] for row in rows) / weighted_denominator if weighted_denominator else 0.0
    def quantile(q):
        return statistics.quantiles(identities, n=100, method="inclusive")[int(q * 100) - 1] if len(identities) > 1 else (identities[0] if identities else 0.0)
    values = [len(rows), weighted_mean, statistics.mean(identities) if identities else 0.0, statistics.median(identities) if identities else 0.0, quantile(.10), quantile(.90), min(identities, default=0.0), max(identities, default=0.0)]
    thresholds = [sum(x >= t for x in identities) for t in (.90, .95, .98, .99)]
    with open(summary_path, "w") as out:
        out.write("sample_id\tprimary_mapped_reads\tidentity_weighted_mean\tidentity_mean\tidentity_median\tidentity_p10\tidentity_p90\tidentity_min\tidentity_max\tidentity_ge_90\tidentity_ge_95\tidentity_ge_98\tidentity_ge_99\n")
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, *values, *thresholds))

if __name__ == "__main__":
    main()
