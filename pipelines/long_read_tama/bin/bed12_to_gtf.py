#!/usr/bin/env python3
"""Convert one or more BED12 model files to sorted exon GTF records."""

import sys


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit("usage: bed12_to_gtf.py OUTPUT.gtf INPUT.bed [...]")
    output, inputs = sys.argv[1], sys.argv[2:]
    records = []
    for path in inputs:
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("track"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    raise SystemExit(f"Invalid BED12 row in {path}: {line.rstrip()}")
                chrom, start, end, name, _, strand = fields[:6]
                block_count = int(fields[9])
                sizes = [int(x) for x in fields[10].rstrip(",").split(",")]
                offsets = [int(x) for x in fields[11].rstrip(",").split(",")]
                if len(sizes) != block_count or len(offsets) != block_count:
                    raise SystemExit(f"Invalid BED12 block list in {path}: {line.rstrip()}")
                transcript_id = name.split(";", 1)[-1].replace('"', "")
                for number, (size, offset) in enumerate(zip(sizes, offsets), 1):
                    exon_start = int(start) + offset + 1
                    exon_end = exon_start + size - 1
                    attrs = f'gene_id "{transcript_id}"; transcript_id "{transcript_id}";'
                    records.append((chrom, exon_start, exon_end, strand, transcript_id, number, attrs))
    records.sort(key=lambda row: (row[0], row[1], row[2], row[4], row[5]))
    with open(output, "w") as handle:
        for chrom, start, end, strand, transcript_id, number, attrs in records:
            handle.write(f"{chrom}\ttmerge\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")


if __name__ == "__main__":
    main()
