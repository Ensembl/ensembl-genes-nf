#!/usr/bin/env python3
"""Prepare a STAR reference with rRNA decoys derived from an Ensembl GFF/GTF."""

import argparse
from pathlib import Path


def fasta_records(path):
    name = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name = line[1:].split()[0]
                sequence = []
            elif name is not None:
                sequence.append(line.strip())
    if name is not None:
        yield name, "".join(sequence)


def attributes(value):
    result = {}
    for field in value.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if " " in field:
            key, item = field.split(None, 1)
            result[key] = item.strip().strip('"')
        elif "=" in field:
            key, item = field.split("=", 1)
            result[key] = item.strip().strip('"')
    return result


def is_rrna(feature, attrs):
    text = " ".join(
        [feature, attrs.get("biotype", ""), attrs.get("gene_biotype", ""),
         attrs.get("gene_type", ""), attrs.get("transcript_biotype", ""),
         attrs.get("transcript_type", "")]
    ).lower()
    return feature.lower() == "rrna" or "rrna" in text or "ribosomal_rna" in text


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-bed", required=True, type=Path)
    parser.add_argument("--output-contigs", required=True, type=Path)
    args = parser.parse_args()

    sequences = dict(fasta_records(args.fasta))
    intervals = []
    with args.gtf.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or not is_rrna(fields[2], attributes(fields[8])):
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            if chrom not in sequences:
                raise SystemExit(f"rRNA annotation contig '{chrom}' is absent from FASTA")
            intervals.append((chrom, start, end))

    decoys = []
    with args.output_fasta.open("w") as output:
        for name, sequence in fasta_records(args.fasta):
            output.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset:offset + 80] + "\n")
        for index, (chrom, start, end) in enumerate(sorted(set(intervals)), 1):
            decoy = f"__rRNA_{args.assembly}_{index}"
            decoys.append((decoy, chrom, start, end, sequences[chrom][start - 1:end]))
            output.write(f">{decoy}\n")
            sequence = decoys[-1][-1]
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset:offset + 80] + "\n")

    with args.output_bed.open("w") as output:
        for decoy, chrom, start, end, _sequence in decoys:
            output.write(f"{decoy}\t0\t{end - start + 1}\t{chrom}:{start}-{end}\n")
    args.output_contigs.write_text("".join(f"{item[0]}\n" for item in decoys))
    if not decoys:
        args.output_bed.write_text("")
        args.output_contigs.write_text("")


if __name__ == "__main__":
    main()
