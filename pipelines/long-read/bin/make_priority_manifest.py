#!/usr/bin/env python3
"""Convert a reviewed priority sheet into a pipeline input manifest.

The discovery sheets deliberately contain raw-subread review candidates.  This
adapter refuses those rows unless they have first been resolved to a supported
PacBio BAM/CCS representation, so a priority sheet cannot be run accidentally.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "id", "cohort_id", "run_id", "technology", "minimap2_preset",
    "isoquant_data_type", "transgram_data_type", "assembly_id", "read_set_id",
    "sample_id", "library_id", "read_representation", "source_uri",
    "selected_artifact_md5", "preprocessing_state", "observation_unit",
    "expected_header_representation", "classification", "proposed_action",
]


def uri(value: str) -> str:
    value = value.strip()
    if value.startswith("ftp."):
        return "ftp://" + value
    return value


def convert(rows: list[dict[str, str]], cohort_id: str, assembly_id: str) -> list[dict[str, str]]:
    output = []
    rejected = []
    for row in rows:
        run_id = row.get("run_accession", "").strip()
        representation = row.get("representation_class", "").strip()
        processing_status = row.get("processing_status", "").strip()
        if processing_status == "REVIEW_ONLY_RAW_SUBREAD_FASTQ" or "RAW_SUBREAD_FASTQ" in representation:
            rejected.append(f"{run_id}: raw-subread FASTQ requires a subreads BAM and PBI")
            continue
        if processing_status not in {"PACBIO_READY_FOR_REVIEW", "ONT_FALLBACK_FASTQ"}:
            rejected.append(f"{run_id}: unsupported processing status {processing_status or 'empty'}")
            continue
        if not run_id or not row.get("url") or not row.get("md5"):
            rejected.append(f"{run_id or '<missing accession>'}: missing URL or MD5")
            continue
        is_pacbio = row.get("platform_priority") in {"1", "2"}
        is_ccs = representation in {"PACBIO_CCS_FASTQ", "PACBIO_CCS_BAM"}
        technology = "pacbio" if is_pacbio else "ont"
        read_representation = "ccs" if is_ccs else "fastq"
        output.append({
            "id": run_id,
            "cohort_id": cohort_id,
            "run_id": run_id,
            "technology": technology,
            "minimap2_preset": "splice:hq" if is_pacbio else "splice",
            "isoquant_data_type": "pacbio_ccs" if is_pacbio else "nanopore",
            "transgram_data_type": "auto",
            "assembly_id": assembly_id,
            "read_set_id": run_id,
            "sample_id": run_id,
            "library_id": run_id,
            "read_representation": read_representation,
            "source_uri": uri(row["url"]),
            "selected_artifact_md5": row["md5"].strip().lower(),
            "preprocessing_state": "ccs_applied" if is_ccs else "processed",
            "observation_unit": "consensus_molecule" if is_pacbio else "ont_read",
            "expected_header_representation": "PACBIO_CCS_ORIGINAL" if is_ccs else "ONT",
            "classification": representation,
            "proposed_action": "ALIGN_PACBIO_CCS" if is_pacbio else "ALIGN_ONT",
        })
    if rejected:
        raise ValueError("priority sheet contains non-runnable rows:\n" + "\n".join(rejected))
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--assembly-id", required=True)
    args = parser.parse_args()
    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    converted = convert(rows, args.cohort_id, args.assembly_id)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(converted)


if __name__ == "__main__":
    main()
