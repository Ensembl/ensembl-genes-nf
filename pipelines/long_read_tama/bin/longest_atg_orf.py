#!/usr/bin/env python3
"""Translate the longest ATG-initiated ORF in each transcript."""

import argparse
from pathlib import Path


STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def fasta_records(path):
    name = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence).upper()
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence).upper()


def longest_orf(sequence):
    candidates = []
    for start in range(len(sequence) - 2):
        if sequence[start:start + 3] != "ATG":
            continue
        stop = None
        for end in range(start + 3, len(sequence) - 2, 3):
            if sequence[end:end + 3] in STOP_CODONS:
                stop = end
                break
        if stop is None:
            end = start + ((len(sequence) - start) // 3) * 3
            status = "PARTIAL_ORF"
        else:
            end = stop + 3
            status = "COMPLETE_ORF"
        coding = sequence[start:end]
        peptide = "".join(CODON_TABLE.get(coding[i:i + 3], "X") for i in range(0, len(coding), 3)).rstrip("*")
        candidates.append((len(peptide), -start, peptide, status, start + 1, end))
    if not candidates:
        return "", "NO_ATG", "", "", ""
    _, _, peptide, status, start, end = max(candidates)
    return peptide, status, start, end, len(peptide)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts")
    parser.add_argument("peptides")
    parser.add_argument("manifest")
    args = parser.parse_args()

    with Path(args.peptides).open("w") as peptide_out, Path(args.manifest).open("w") as manifest_out:
        manifest_out.write("model_id\ttranscript_id\tpeptide_id\torf_status\torf_start\torf_end\tpeptide_length\n")
        for model_id, sequence in fasta_records(Path(args.transcripts)):
            peptide, status, start, end, length = longest_orf(sequence)
            peptide_id = f"{model_id}.orf1" if peptide else ""
            manifest_out.write(f"{model_id}\t{model_id}\t{peptide_id}\t{status}\t{start}\t{end}\t{length}\n")
            if peptide:
                peptide_out.write(f">{peptide_id} model_id={model_id}\n{peptide}\n")


if __name__ == "__main__":
    main()
