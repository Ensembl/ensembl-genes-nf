#!/usr/bin/env python3
"""Convert GenBlast protein-alignment GTFs to Funannotate GFF3."""

import argparse
import re
from pathlib import Path

TRANSCRIPT_ID = re.compile(r'transcript_id\s+"([^"]+)"')


def convert(inputs, output):
    with output.open("w") as destination:
        for input_path in inputs:
            tag = re.sub(r"^protein_alignment_\d+_", "", input_path.stem)
            transcript_number = 0
            amino_acids = 0
            alignment_id = None
            protein_id = None
            score = "."

            with input_path.open() as source:
                for line_number, line in enumerate(source, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        raise ValueError(f"{input_path}:{line_number}: expected 9 GTF fields")

                    feature = fields[2]
                    if feature == "transcript":
                        transcript_number += 1
                        score = fields[5]
                        match = TRANSCRIPT_ID.search(fields[8])
                        if not match:
                            raise ValueError(
                                f"{input_path}:{line_number}: transcript_id attribute is missing"
                            )
                        protein_id = match.group(1)
                        alignment_id = f"{tag}_{transcript_number}_{protein_id}"
                        amino_acids = 0
                    elif feature == "exon":
                        if alignment_id is None or protein_id is None:
                            raise ValueError(f"{input_path}:{line_number}: exon occurs before a transcript")
                        start, end = int(fields[3]), int(fields[4])
                        length = end - start + 1
                        protein_start = amino_acids + 1
                        amino_acids += (length + 2) // 3
                        attributes = (
                            f"ID={alignment_id};Target={protein_id} " f"{protein_start} {amino_acids}"
                        )
                        destination.write(
                            "\t".join(
                                [
                                    fields[0],
                                    "genblast",
                                    "nucleotide_to_protein_match",
                                    fields[3],
                                    fields[4],
                                    score,
                                    fields[6],
                                    ".",
                                    attributes,
                                ]
                            )
                            + "\n"
                        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    convert(args.inputs, args.output)


if __name__ == "__main__":
    main()
