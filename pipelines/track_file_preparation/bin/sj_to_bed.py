#!/usr/bin/env python3
"""Convert STAR SJ.out.tab records to the track-prep splice BED schema."""

from __future__ import annotations

import sys
from pathlib import Path


MOTIFS = {
    0: ".",
    1: "GT-AG",
    2: "CT-AC",
    3: "GC-AG",
    4: "CT-GC",
    5: "AT-AC",
    6: "GT-AT",
    7: "GC-AT",
}
STRANDS = {0: ".", 1: "+", 2: "-"}


def convert(source: Path, destination: Path) -> None:
    records = 0
    with source.open() as src, destination.open("w") as out:
        for line_number, line in enumerate(src, start=1):
            fields = line.rstrip("\n").split("\t")
            if line_number == 1 and fields[:9] == [
                "chrom", "start", "end", "strand", "motif", "annotated",
                "uniqReads", "multiReads", "maxOverhang",
            ]:
                continue
            if len(fields) != 9:
                raise ValueError(f"{source}:{line_number}: expected 9 STAR columns, got {len(fields)}")
            chrom = fields[0]
            try:
                start, end, strand_code, motif_code, unique_reads, multi_reads = (
                    int(fields[index]) for index in (1, 2, 3, 4, 6, 7)
                )
            except ValueError as error:
                raise ValueError(f"{source}:{line_number}: non-integer STAR field") from error
            if not chrom or start < 1 or end < start:
                raise ValueError(f"{source}:{line_number}: invalid chromosome or junction coordinates")
            if strand_code not in STRANDS or motif_code not in MOTIFS:
                raise ValueError(f"{source}:{line_number}: unsupported strand or STAR motif code")
            if unique_reads < 0 or multi_reads < 0:
                raise ValueError(f"{source}:{line_number}: read counts must be non-negative")
            strand = STRANDS[strand_code]
            name = f"{chrom}:{start}-{end}:{strand}"
            score = min(unique_reads, 1000)
            out.write("\t".join(map(str, (
                chrom, start - 1, end, name, score, strand,
                unique_reads, multi_reads, MOTIFS[motif_code],
            ))) + "\n")
            records += 1
    if records == 0:
        raise ValueError(f"{source}: no STAR splice-junction records found")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: sj_to_bed.py STAR_SJ.out.tab output.bed")
    convert(Path(sys.argv[1]), Path(sys.argv[2]))


if __name__ == "__main__":
    main()
