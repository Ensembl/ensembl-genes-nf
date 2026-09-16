#!/usr/bin/env python3
"""Create tiny, deterministic local inputs for the long-read TAMA examples."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--url", required=True, help="HTTP URL prefix served by run_local_example.sh")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    (args.outdir / "genome.fa").write_text(
        ">chr1\n" + "ACGT" * 75 + "\n"
        ">chr2\n" + "TGCA" * 75 + "\n"
        ">chr3\n" + "GATC" * 75 + "\n"
    )
    fastq = args.outdir / "local_ont.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        for number, sequence in enumerate(("ACGT" * 20, "TGCA" * 20), 1):
            read_id = f"123e4567-e89b-12d3-a456-42661417400{number}"
            handle.write(f"@{read_id} runid=local-example ch=1\n{sequence}\n+\n" + "I" * len(sequence) + "\n")

    md5 = hashlib.md5(fastq.read_bytes()).hexdigest()
    url = args.url.rstrip("/") + "/" + fastq.name
    (args.outdir / "metadata.json").write_text(json.dumps({"SRR900001": {
        "instrument_platform": "ONT", "instrument_model": "MinION",
        "sra_platform": "ONT", "fastq_ftp": url, "fastq_md5": md5,
        "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC"
    }}, indent=2) + "\n")
    with (args.outdir / "candidate_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_accession", "tissue", "description", "url", "md5", "platform", "filename"])
        writer.writeheader()
        writer.writerow({"run_accession": "SRR900001", "tissue": "local-test", "description": "two-read ONT smoke fixture",
                         "url": url, "md5": md5, "platform": "ONT", "filename": fastq.name})

    # This is deliberately a reviewed input example. The pipeline validates
    # it, then exercises acquisition, alignment, collapse, merge, and QC in
    # stub mode without requiring a site-specific TAMA installation.
    fields = ["run_accession", "tissue", "description", "classification", "proposed_action",
              "selected_artifact_uri", "selected_artifact_md5", "selected_artifact_basename",
              "minimap2_preset", "expected_header_representation", "classification_report_sha256",
              "reviewer", "reviewed_at", "status", "review_decision"]
    with (args.outdir / "approved_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({"run_accession": "SRR900001", "tissue": "local-test", "description": "two-read ONT smoke fixture",
                         "classification": "ONT_FASTQ", "proposed_action": "ALIGN_ONT", "selected_artifact_uri": url,
                         "selected_artifact_md5": md5, "selected_artifact_basename": fastq.name, "minimap2_preset": "splice",
                         "expected_header_representation": "ONT", "classification_report_sha256": "0" * 64,
                         "reviewer": "local-example", "reviewed_at": "2026-09-16T00:00:00Z", "status": "APPROVED", "review_decision": "APPROVE"})


if __name__ == "__main__":
    main()
