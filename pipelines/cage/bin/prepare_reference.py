#!/usr/bin/env python3
"""Prepare a STAR reference with rRNA decoys derived from an Ensembl GFF/GTF."""

import argparse
from pathlib import Path


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

    intervals = {}
    with args.gtf.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or not is_rrna(fields[2], attributes(fields[8])):
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            intervals.setdefault(chrom, set()).add((start, end))

    decoys = []
    sequence_parts = {}
    with args.output_fasta.open("w") as output:
        current = None
        offset = 0
        with args.fasta.open() as fasta:
            for raw_line in fasta:
                line = raw_line.rstrip("\n")
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    offset = 0
                    sequence_parts[current] = {}
                    output.write(raw_line)
                    continue
                if current is None:
                    continue
                sequence = line.strip()
                output.write(raw_line)
                for start, end in intervals.get(current, ()):
                    left = max(start - 1 - offset, 0)
                    right = min(end - offset, len(sequence))
                    if left < right:
                        sequence_parts[current].setdefault((start, end), []).append(sequence[left:right])
                offset += len(sequence)

        missing = sorted(set(intervals) - set(sequence_parts))
        if missing:
            raise SystemExit(f"rRNA annotation contig(s) absent from FASTA: {', '.join(missing[:5])}")

        index = 0
        for chrom in sorted(intervals):
            for start, end in sorted(intervals[chrom]):
                index += 1
                sequence = "".join(sequence_parts[chrom][(start, end)])
                decoy = f"__rRNA_{args.assembly}_{index}"
                decoys.append((decoy, chrom, start, end, sequence))
                output.write(f">{decoy}\n")
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
