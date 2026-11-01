#!/usr/bin/env python3
"""Convert exon GTF records into deterministic, strand-preserving BED12."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ATTR = re.compile(r'(\w+)\s+"([^"]*)"')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--bed12", required=True, type=Path)
    parser.add_argument("--backend", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    transcripts = defaultdict(list)
    with args.gtf.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(ATTR.findall(fields[8]))
            tid = attrs.get("transcript_id") or attrs.get("gene_id")
            if not tid:
                raise SystemExit("exon record is missing transcript_id")
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise SystemExit(f"invalid GTF interval for {tid}")
            transcripts[(fields[0], tid, fields[6] if fields[6] in "+-" else ".")].append((start, end))
    if not transcripts:
        if any(line.startswith("# EMPTY") for line in args.gtf.read_text().splitlines()):
            args.bed12.touch()
            write_auxiliary_outputs(args, 0)
            return
        raise SystemExit("GTF contains no exon records")
    with args.bed12.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for index, ((chrom, tid, strand), exons) in enumerate(sorted(transcripts.items()), 1):
            exons = sorted(set(exons))
            chrom_start, chrom_end = exons[0][0], exons[-1][1]
            sizes = ",".join(str(end - start) for start, end in exons) + ","
            starts = ",".join(str(start - chrom_start) for start, _end in exons) + ","
            name = f"{args.backend}.{args.scope}.model_{index:06d}.{tid}"
            writer.writerow([chrom, chrom_start, chrom_end, name, 0, strand, chrom_start, chrom_end, 0, len(exons), sizes, starts])
    write_auxiliary_outputs(args, len(transcripts))


def write_auxiliary_outputs(args: argparse.Namespace, model_count: int) -> None:
    if args.report:
        checksum = hashlib.sha256(args.bed12.read_bytes()).hexdigest()
        args.report.write_text(f"metric\tvalue\nmodels\t{model_count}\nerrors\t0\nsha256\t{checksum}\n")
    if args.manifest:
        args.manifest.write_text(json.dumps({"backend": args.backend, "scope": args.scope, "source": str(args.gtf), "models": model_count}, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
