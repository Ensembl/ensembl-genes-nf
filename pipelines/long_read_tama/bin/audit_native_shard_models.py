#!/usr/bin/env python3
"""Count native per-contig transcript models and retain only mergeable shards."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


TRANSCRIPT_ID = re.compile(r'(?:^|;)\s*transcript_id\s+"([^"]+)"')
GFF_ID = re.compile(r'(?:^|;)\s*ID=([^;]+)')


def count_gtf_models(path: Path) -> int:
    transcript_ids: set[str] = set()
    transcript_records = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 GTF columns, found {len(fields)}")
            feature = fields[2].lower()
            # Gene and exon rows are not transcript models and may not carry
            # transcript_id attributes.
            if feature != "transcript":
                continue
            attributes = fields[8]
            match = TRANSCRIPT_ID.search(attributes)
            if match is None:
                match = GFF_ID.search(attributes)
            if match is None:
                raise ValueError(f"{path}:{line_number}: transcript feature has no transcript identifier")
            transcript_records += 1
            transcript_ids.add(match.group(1))
    return len(transcript_ids) if transcript_ids else transcript_records


def count_bed12_models(path: Path) -> int:
    model_records = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{line_number}: expected at least 12 BED12 columns, found {len(fields)}")
            model_records += 1
    return model_records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model", type=Path)
    parser.add_argument("native_format", choices=("gtf", "bed12"))
    parser.add_argument("backend")
    parser.add_argument("accession")
    parser.add_argument("shard")
    parser.add_argument("mergeable_model", type=Path)
    parser.add_argument("counts", type=Path)
    args = parser.parse_args()

    if not args.model.is_file():
        raise SystemExit(f"native model file is missing: {args.model}")
    if args.native_format == "gtf":
        model_count = count_gtf_models(args.model)
    else:
        model_count = count_bed12_models(args.model)

    status = "MERGEABLE" if model_count else "EMPTY"
    if model_count:
        shutil.copyfile(args.model, args.mergeable_model)

    with args.counts.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("backend", "run_accession", "shard", "native_format", "model_count", "merge_status"))
        writer.writerow((args.backend, args.accession, args.shard, args.native_format, model_count, status))


if __name__ == "__main__":
    main()
