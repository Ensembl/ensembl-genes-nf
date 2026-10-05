#!/usr/bin/env python3
"""Audit stable IsoQuant products and write a complete/failed status manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


def non_header_records(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#")) - 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("products")
    parser.add_argument("output")
    parser.add_argument("status")
    args = parser.parse_args()
    root = Path(args.products)
    provenance = root / "isoquant_input_manifest.tsv"
    required = [root / "transcript_models.gtf", root / "transcript_model_reads.tsv.gz", root / "read_info.tsv.gz"]
    missing = [str(path.name) for path in required if not path.is_file() or path.stat().st_size == 0]
    model_count = 0
    read_info_count = 0
    read_to_model_count = 0
    if not missing:
        model_count = sum(1 for line in required[0].read_text().splitlines() if line.strip() and not line.startswith("#") and "\ttranscript\t" in line)
        read_info_count = max(non_header_records(required[2]), 0)
        read_to_model_count = max(non_header_records(required[1]), 0)
    status = "COMPLETE" if not missing and model_count > 0 and read_info_count > 0 and read_to_model_count > 0 else "EMPTY_RESULT"
    if missing:
        status = "FAILED_PARTIAL"
    rows = {}
    if provenance.exists():
        with provenance.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows[row["key"]] = row["value"]
    read_info_source = "unknown"
    product_manifest = root / "product_manifest.tsv"
    if product_manifest.exists():
        with product_manifest.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["stable_name"] == "read_info.tsv.gz":
                    read_info_source = row.get("source_type", "unknown")
    output_sha = sha256(required[0]) if required[0].exists() else "NONE"
    fields = [
        ("backend", "isoquant"),
        ("scope", rows.get("scope", "UNKNOWN")),
        ("accession_or_cohort", rows.get("accession_or_cohort", "UNKNOWN")),
        ("isoquant_version", rows.get("isoquant_version", "UNKNOWN")),
        ("data_type", rows.get("data_type", "UNKNOWN")),
        ("mode", rows.get("mode", "UNKNOWN")),
        ("reference_fasta_sha256", rows.get("reference_fasta_sha256", "UNKNOWN")),
        ("genedb_sha256_or_NONE", rows.get("genedb_sha256_or_NONE", "NONE")),
        ("input_bam_names", rows.get("input_bam_names", "UNKNOWN")),
        ("input_bam_checksums", rows.get("input_bam_checksums", "UNKNOWN")),
        ("transcript_model_count", str(model_count)),
        ("read_info_count", str(read_info_count)),
        ("read_to_model_record_count", str(read_to_model_count)),
        ("read_info_source", read_info_source),
        ("output_file", required[0].name),
        ("output_sha256", output_sha),
        ("complete", status == "COMPLETE"),
        ("status", status),
        ("missing_products", ",".join(missing) or "NONE"),
    ]
    with Path(args.output).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(fields)
    Path(args.status).write_text(f"status\tcomplete\tmodels\tread_info\tread_to_model\n{status}\t{status == 'COMPLETE'}\t{model_count}\t{read_info_count}\t{read_to_model_count}\n")
    if status != "COMPLETE":
        raise SystemExit(f"IsoQuant output is {status}; products={root}")


if __name__ == "__main__":
    main()
