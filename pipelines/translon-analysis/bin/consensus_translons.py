#!/usr/bin/env python3
"""Build a deterministic cross-caller candidate-translon contract."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


CODONS = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def read_fasta(path):
    sequences, name, chunks = {}, None, []
    for line in path.open():
        line = line.strip()
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        elif line:
            chunks.append(line)
    if name:
        sequences[name] = "".join(chunks).upper()
    return sequences


def peptide(row, genome):
    sequence = genome.get(row["chrom"])
    if not sequence:
        return ""
    dna = sequence[row["start"]:row["end"]]
    if row["strand"] == "-":
        dna = dna.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    try:
        frame = int(row["frame"])
    except ValueError:
        frame = 0
    dna = dna[frame:]
    return "".join(CODONS.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


FIELDS = ["interval_id", "chrom", "start", "end", "strand", "frame",
          "caller_count", "callers", "status"]


def read_rows(paths):
    for path in sorted(paths):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if not row.get("chrom") or not row.get("start") or not row.get("end"):
                    continue
                try:
                    start, end = int(row["start"]), int(row["end"])
                except ValueError:
                    continue
                if end <= start:
                    continue
                tool = row.get("tool") or path.stem.split(".")[0]
                yield {
                    "chrom": row["chrom"], "start": start, "end": end,
                    "strand": row.get("strand") or ".",
                    "frame": row.get("frame") or ".", "tool": tool,
                    "sample_id": row.get("sample_id", "unknown"),
                    "orf_id": row.get("orf_id", ""),
                    "score": row.get("score", ""),
                }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--min-caller-agreement", type=int, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--bed12", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument("--verdicts", type=Path, required=True)
    args = parser.parse_args()
    genome = read_fasta(args.genome_fasta)

    groups = defaultdict(list)
    for row in read_rows(args.input_dir.glob("*.tsv")):
        key = (row["chrom"], row["start"], row["end"], row["strand"], row["frame"])
        groups[key].append(row)

    consensus = []
    for n, (key, rows) in enumerate(sorted(groups.items()), 1):
        chrom, start, end, strand, frame = key
        callers = sorted({r["tool"] for r in rows})
        interval_id = f"translon_{n:07d}"
        consensus.append({
            "interval_id": interval_id, "chrom": chrom, "start": start,
            "end": end, "strand": strand, "frame": frame,
            "caller_count": len(callers), "callers": ",".join(callers),
            "status": "trusted" if len(callers) >= args.min_caller_agreement else "candidate",
        })

    with args.candidates.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader(); writer.writerows(consensus)

    with args.intervals.open("w", newline="") as handle:
        interval_fields = ["interval_id", "chrom", "start", "end", "strand", "frame"]
        writer = csv.DictWriter(handle, fieldnames=interval_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: row[field] for field in interval_fields}
                         for row in consensus if row["status"] == "trusted")

    with args.verdicts.open("w", newline="") as handle:
        fields = ["interval_id", "frame", "mechanism_class", "confidence", "condition_state", "peptide_sequence"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in consensus:
            if row["status"] == "trusted":
                writer.writerow({"interval_id": row["interval_id"], "frame": row["frame"],
                                 "mechanism_class": "caller_consensus",
                                 "confidence": row["caller_count"],
                                 "condition_state": "condition-unresolved",
                                 "peptide_sequence": peptide(row, genome)})

    with args.bed12.open("w") as handle:
        for row in consensus:
            block = row["end"] - row["start"]
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['interval_id']}\t{row['caller_count']}\t{row['strand']}\t{row['start']}\t{row['end']}\t0\t1\t{block},\t0,\n")


if __name__ == "__main__":
    main()
