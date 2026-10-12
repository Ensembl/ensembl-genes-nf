#!/usr/bin/env python3
"""Read GFF3/GTF once and write a normalized, tab-separated feature file."""

from __future__ import annotations

import argparse
import csv
import gzip
from urllib.parse import unquote

from qc_utils import (
    BIOTYPE_KEYS, NORMALIZED_FIELDS, first, is_canonical, open_text,
    parse_gff3_attributes, parse_gtf_attributes, split_annotation_line, values,
)


TRANSCRIPT_TYPES = {"transcript", "mrna", "rrna", "trna", "ncrna", "lncrna", "snrna", "snorna"}


def detect_format(path: str) -> str:
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.startswith("##gff-version"):
                return "gff3"
            if line.startswith("#"):
                continue
            columns = split_annotation_line(line)
            if len(columns) >= 9:
                return "gtf" if '"' in columns[8] and "=" not in columns[8] else "gff3"
    return "gtf" if str(path).lower().endswith((".gtf", ".gtf.gz")) else "gff3"


def normalize(path: str, annotation_format: str):
    annotation_format = detect_format(path) if annotation_format == "auto" else annotation_format
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line or line.startswith("#"):
                continue
            columns = split_annotation_line(line)
            if len(columns) < 9:
                continue
            try:
                attrs = (
                    parse_gtf_attributes(columns[8])
                    if annotation_format == "gtf"
                    else parse_gff3_attributes(columns[8])
                )
                feature_type = columns[2]
                feature_id = first(attrs, ("gene_id", "transcript_id", "ID", "Name"))
                parent_ids = values(attrs, ("Parent",))
                transcript_ids = values(attrs, ("transcript_id",))
                if not transcript_ids and feature_type.lower() in {"exon", "cds"}:
                    transcript_ids = parent_ids
                gene_id = first(attrs, ("gene_id", "gene"))
                if not gene_id and feature_type.lower() == "gene":
                    gene_id = feature_id
                biotype = first(attrs, BIOTYPE_KEYS)
                yield {
                    "seqid": unquote(columns[0]), "source": columns[1],
                    "feature_type": feature_type, "start": int(columns[3]) - 1,
                    "end": int(columns[4]), "score": columns[5], "strand": columns[6],
                    "phase": columns[7], "feature_id": feature_id,
                    "gene_id": gene_id, "parent_ids": ",".join(parent_ids),
                    "transcript_ids": ",".join(transcript_ids), "biotype": biotype,
                    "is_canonical": "1" if is_canonical(attrs) else "0",
                }
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Invalid annotation row {line_number}: {exc}") from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotation")
    parser.add_argument("--format", choices=("auto", "gff3", "gtf"), default="auto")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open_text(args.output, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalize(args.annotation, args.format))


if __name__ == "__main__":
    main()


