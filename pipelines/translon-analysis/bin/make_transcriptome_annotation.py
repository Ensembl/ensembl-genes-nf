#!/usr/bin/env python3
"""Project a genome-space GTF/FASTA into transcript-space coordinates."""
import argparse
import collections
import re
from pathlib import Path

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def attrs(text):
    out = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', text):
        out[key] = value
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out-gtf", required=True, type=Path)
    parser.add_argument("--out-fasta", required=True, type=Path)
    args = parser.parse_args()

    records = collections.defaultdict(list)
    with args.gtf.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            record = attrs(fields[8])
            tx = record.get("transcript_id")
            if tx:
                records[tx].append((fields, record))

    fasta = pysam.FastaFile(str(args.fasta))
    args.out_gtf.parent.mkdir(parents=True, exist_ok=True)
    with args.out_gtf.open("w") as gtf, args.out_fasta.open("w") as out_fa:
        for tx, entries in records.items():
            exons = [x for x in entries if x[0][2].lower() == "exon"]
            if not exons:
                continue
            strand = exons[0][0][6]
            exons.sort(key=lambda x: int(x[0][3]), reverse=strand == "-")
            tx_len = sum(int(x[0][4]) - int(x[0][3]) + 1 for x in exons)
            sequence = "".join(
                fasta.fetch(x[0][0], int(x[0][3]) - 1, int(x[0][4])).upper()
                for x in exons
            )
            if strand == "-":
                sequence = sequence.translate(COMPLEMENT)[::-1]
            out_fa.write(f">{tx}\n")
            for i in range(0, len(sequence), 80):
                out_fa.write(sequence[i : i + 80] + "\n")

            offsets = []
            cursor = 1
            for fields, _ in exons:
                start, end = int(fields[3]), int(fields[4])
                offsets.append((fields[0], start, end, cursor))
                cursor += end - start + 1
            gene_id = next((record.get("gene_id", tx) for _, record in entries), tx)
            gene_name = next((record.get("gene_name", gene_id) for _, record in entries), gene_id)
            gtf.write(
                "\t".join([tx, "translon", "gene", "1", str(tx_len), ".", "+", ".",
                           f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{tx}";'])
                + "\n"
            )
            for fields, record in entries:
                feature_start, feature_end = int(fields[3]), int(fields[4])
                for chrom, exon_start, exon_end, tx_offset in offsets:
                    start = max(feature_start, exon_start)
                    end = min(feature_end, exon_end)
                    if fields[0] != chrom or start > end:
                        continue
                    if strand == "+":
                        rel_start = tx_offset + start - exon_start
                        rel_end = tx_offset + end - exon_start
                    else:
                        rel_start = tx_offset + exon_end - end
                        rel_end = tx_offset + exon_end - start
                    projected = list(fields)
                    projected[0] = tx
                    projected[3], projected[4] = str(rel_start), str(rel_end)
                    projected[6] = "+"
                    projected[8] = "; ".join(f'{k} "{v}"' for k, v in record.items()) + ";"
                    gtf.write("\t".join(projected) + "\n")


if __name__ == "__main__":
    main()
