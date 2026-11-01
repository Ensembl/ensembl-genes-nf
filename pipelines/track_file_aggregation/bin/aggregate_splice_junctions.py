#!/usr/bin/env python3
"""Sum splice-junction support while preserving the shared nine-field schema."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP


def weighted_count(value: int, divisor: Decimal) -> int:
    return int((Decimal(value) / divisor).quantize(Decimal("1"), rounding=ROUND_HALF_UP))


def aggregate(inputs: list[Path], output: Path, weight_divisor: Decimal = Decimal("1")) -> None:
    if weight_divisor <= 0:
        raise ValueError("weight divisor must be greater than zero")
    records: dict[tuple[str, int, int, str], dict[str, object]] = {}
    motifs: defaultdict[tuple[str, int, int, str], set[str]] = defaultdict(set)
    for path in inputs:
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{line_number}: expected at least nine BED fields")
            chrom, start, end, _name, _score, strand = fields[:6]
            unique_reads, multi_reads = fields[6:8]
            motif = fields[8]
            key = (chrom, int(start), int(end), strand)
            motifs[key].add(motif)
            records.setdefault(key, {"unique_reads": 0, "multi_reads": 0})
            records[key]["unique_reads"] = int(records[key]["unique_reads"]) + int(unique_reads)
            records[key]["multi_reads"] = int(records[key]["multi_reads"]) + int(multi_reads)

    output_rows: list[tuple[str, int, int, str, int, str, int, int, str]] = []
    for key in sorted(records, key=lambda item: (item[0], item[1], item[2], item[3])):
        chrom, start, end, strand = key
        known = {motif for motif in motifs[key] if motif != "."}
        if len(known) > 1:
            raise ValueError(f"conflicting splice motifs for {chrom}:{start}-{end}:{strand}: {sorted(known)}")
        motif = next(iter(known), ".")
        unique_reads = weighted_count(int(records[key]["unique_reads"]), weight_divisor)
        multi_reads = weighted_count(int(records[key]["multi_reads"]), weight_divisor)
        name = f"{chrom}:{start + 1}-{end}:{strand}"
        output_rows.append((chrom, start, end, name, min(unique_reads, 1000), strand, unique_reads, multi_reads, motif))

    with output.open("w", encoding="utf-8") as handle:
        for row in output_rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--weight-divisor", type=Decimal, default=Decimal("1"))
    args = parser.parse_args()
    aggregate(args.inputs, args.output, args.weight_divisor)


if __name__ == "__main__":
    main()
