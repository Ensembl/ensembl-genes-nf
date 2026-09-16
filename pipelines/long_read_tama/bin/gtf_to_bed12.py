#!/usr/bin/env python3
"""Convert exon-only GTF models into BED12."""

import re
import sys
from collections import defaultdict


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: gtf_to_bed12.py INPUT.gtf OUTPUT.bed COHORT")
    source, output, cohort = sys.argv[1:]
    models = defaultdict(list)
    with open(source) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = re.search(r'transcript_id\s+"([^"]+)"', fields[8])
            if not match:
                raise SystemExit(f"Missing transcript_id in GTF row: {line.rstrip()}")
            models[match.group(1)].append(fields)
    rows = []
    for transcript_id, exons in models.items():
        exons.sort(key=lambda row: int(row[3]))
        chrom = exons[0][0]
        strand = exons[0][6]
        if any(row[0] != chrom or row[6] != strand for row in exons):
            raise SystemExit(f"Transcript {transcript_id} spans chromosomes or strands")
        chrom_start = min(int(row[3]) for row in exons) - 1
        chrom_end = max(int(row[4]) for row in exons)
        sizes = [int(row[4]) - int(row[3]) + 1 for row in exons]
        offsets = [int(row[3]) - 1 - chrom_start for row in exons]
        name = f"{cohort};{transcript_id}"
        rows.append((chrom, chrom_start, chrom_end, name, strand, len(exons), sizes, offsets))
    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    with open(output, "w") as handle:
        for chrom, start, end, name, strand, count, sizes, offsets in rows:
            handle.write("\t".join([
                chrom, str(start), str(end), name, "0", strand, str(start), str(end),
                "0", str(count), ",".join(map(str, sizes)) + ",", ",".join(map(str, offsets)) + ",",
            ]) + "\n")


if __name__ == "__main__":
    main()
