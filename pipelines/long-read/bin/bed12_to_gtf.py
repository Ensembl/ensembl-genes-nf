#!/usr/bin/env python3
"""Convert TAMA BED12 records to deterministic exon-level GTF."""

import argparse
from pathlib import Path
from typing import List, Tuple


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--source", default="pipeline")
    args = parser.parse_args()

    records: List[Tuple[str, int, str, List[Tuple[int, int]], str]] = []
    with args.bed12.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                raise SystemExit(f"BED12 line {line_number} has {len(fields)} columns")
            chrom, start_text, end_text, name, _score, strand = fields[:6]
            start, end = int(start_text), int(end_text)
            block_count = int(fields[9])
            sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
            offsets = [int(value) for value in fields[11].rstrip(",").split(",") if value]
            if block_count != len(sizes) or block_count != len(offsets) or end <= start:
                raise SystemExit(f"invalid BED12 block lists on line {line_number}")
            if strand not in {"+", "-", "."}:
                raise SystemExit(f"invalid BED12 strand on line {line_number}: {strand}")
            exons = []
            for size, offset in zip(sizes, offsets):
                exon_start = start + offset
                exon_end = exon_start + size
                if size <= 0 or exon_start < start or exon_end > end:
                    raise SystemExit(f"BED12 block outside transcript bounds on line {line_number}")
                exons.append((exon_start, exon_end))
            records.append((chrom, start, strand, exons, name))

    records.sort(key=lambda record: (record[0], record[1], record[4], record[2]))
    with args.gtf.open("w") as handle:
        for _chrom, _start, _strand, exons, name in records:
            transcript_id = name
            attributes = f'gene_id "{transcript_id}"; transcript_id "{transcript_id}";'
            for exon_number, (exon_start, exon_end) in enumerate(exons, 1):
                fields = [
                    _chrom,
                    args.source,
                    "exon",
                    exon_start + 1,
                    exon_end,
                    ".",
                    _strand,
                    ".",
                    f'{attributes} exon_number "{exon_number}";',
                ]
                handle.write("\t".join(map(str, fields)) + "\n")

    if not records:
        raise SystemExit("TAMA BED12 contains no transcript records")


if __name__ == "__main__":
    main()
