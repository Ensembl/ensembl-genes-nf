#!/usr/bin/env python3
"""Promote safe selected runs and audit excluded candidate rows."""

import csv
import hashlib
import sys
from datetime import datetime, timezone

FIELDS = ["run_accession", "tissue", "description", "classification", "proposed_action",
          "selected_artifact_uri", "selected_artifact_md5", "selected_artifact_basename",
          "minimap2_preset", "expected_header_representation", "classification_report_sha256",
          "reviewer", "reviewed_at", "status", "review_decision"]
ALLOWED = {"PACBIO_CCS_FASTQ", "PACBIO_PROCESSED_FASTQ", "PACBIO_CCS_BAM", "PACBIO_SUBREAD_BAM", "ONT_FASTQ"}


def main():
    report_dir, output, audit_output = sys.argv[1:]
    report = f"{report_dir}/run_classification.tsv"
    report_hash = hashlib.sha256(open(report, "rb").read()).hexdigest()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    rows = []
    audit_rows = []
    with open(report, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("run_accession", "")
            classification = row.get("classification", "")
            status = row.get("status", "")
            if status != "READY_FOR_REVIEW" or classification not in ALLOWED:
                audit_rows.append({"run_accession": accession, "classification": classification,
                                   "status": status, "decision": "QUARANTINE",
                                   "detail": row.get("reason_codes", "not eligible for automatic processing")})
                continue
            expected = "ONT" if classification == "ONT_FASTQ" else (
                "PACBIO_CCS_ORIGINAL" if "NCBI_ORIGINAL_CCS_FASTQ" in row.get("reason_codes", "") else
                "PACBIO_CCS" if classification in {"PACBIO_CCS_FASTQ", "PACBIO_CCS_BAM", "PACBIO_SUBREAD_BAM"} else "PACBIO_PROCESSED"
            )
            rows.append({
                "run_accession": row["run_accession"], "tissue": row.get("tissue", "unknown"),
                "description": row.get("description", "unknown"), "classification": classification,
                "proposed_action": row["proposed_action"], "selected_artifact_uri": row["selected_artifact_uri"],
                "selected_artifact_md5": row["selected_artifact_md5"], "selected_artifact_basename": row["selected_artifact_basename"],
                "minimap2_preset": "splice" if classification == "ONT_FASTQ" else "splice:hq",
                "expected_header_representation": expected, "classification_report_sha256": report_hash,
                "reviewer": "automatic_selector", "reviewed_at": reviewed_at,
                "status": "APPROVED", "review_decision": "APPROVE",
            })
            audit_rows.append({"run_accession": accession, "classification": classification,
                               "status": "APPROVED", "decision": "PROCESS",
                               "detail": "READY_FOR_REVIEW and compatible classification"})
    if not rows:
        raise SystemExit("automatic selection produced no safely runnable runs")
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with open(audit_output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_accession", "classification", "status", "decision", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)


if __name__ == "__main__":
    main()
