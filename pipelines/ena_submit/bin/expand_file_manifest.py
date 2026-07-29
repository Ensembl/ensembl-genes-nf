#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from pathlib import Path


OUT_COLUMNS = [
    "analysis_id",
    "project_alias",
    "assembly",
    "release",
    "study",
    "umbrella_study",
    "analysis_alias",
    "title",
    "description",
    "assembly_accession",
    "reference_fasta",
    "assembly_report",
    "last_geneset_update",
    "partial_release_label",
    "species",
    "taxon_id",
    "ref_seqs",
    "analysis_links",
    "analysis_attributes",
    "analysis_type",
    "omit_run_refs_in_test",
    "file_path",
    "file_type",
    "remote_name",
    "run_accession",
    "sample_accession",
    "experiment_accession",
]


def read_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def slug(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return re.sub(r"_+", "_", value).strip("_")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate and normalize one annotation files.tsv for Nextflow.")
    parser.add_argument("--files-tsv", required=True)
    parser.add_argument("--analysis-tsv", required=True)
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--project-alias", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--default-file-type", default="")
    parser.add_argument("--out", default="expanded_files.tsv")
    args = parser.parse_args()

    analysis_rows = read_rows(Path(args.analysis_tsv))
    if len(analysis_rows) != 1:
        raise RuntimeError(f"analysis-tsv must contain exactly one row: {args.analysis_tsv}")
    analysis = analysis_rows[0]

    rows = read_rows(Path(args.files_tsv))
    if not rows:
        raise RuntimeError(f"files_tsv must contain at least one file row: {args.files_tsv}")

    out_rows = []
    for idx, row in enumerate(rows, start=2):
        file_path = (row.get("file_path") or "").strip()
        if not file_path:
            raise RuntimeError(f"{args.files_tsv}:{idx}: missing file_path")
        if not os.path.exists(file_path):
            raise RuntimeError(f"{args.files_tsv}:{idx}: missing alignment file: {file_path}")

        file_type = (row.get("file_type") or args.default_file_type or "").strip().lower()
        if file_type not in {"bam", "cram"}:
            raise RuntimeError(f"{args.files_tsv}:{idx}: file_type must be bam or cram")

        remote_name = (row.get("remote_name") or "").strip() or os.path.basename(file_path)
        run_accession = (row.get("run_accession") or row.get("run_accessions") or "").strip()
        analysis_suffix = slug(run_accession or Path(remote_name).stem)
        file_analysis_alias = slug(f"{analysis.get('analysis_alias', args.analysis_id)}_{analysis_suffix}")
        out_rows.append(
            {
                # ENA accepts one BAM or CRAM per ANALYSIS. The input manifest
                # remains annotation-level, but each expanded file gets its
                # own stable analysis identity.
                "analysis_id": file_analysis_alias,
                "project_alias": args.project_alias,
                "assembly": args.assembly,
                "release": args.release,
                "study": analysis.get("study", ""),
                "umbrella_study": analysis.get("umbrella_study", ""),
                "analysis_alias": file_analysis_alias,
                "title": f"{analysis.get('title', '')} ({run_accession or analysis_suffix})",
                "description": analysis.get("description", ""),
                "assembly_accession": analysis.get("assembly_accession", ""),
                "reference_fasta": analysis.get("reference_fasta", ""),
                "assembly_report": analysis.get("assembly_report", ""),
                "last_geneset_update": analysis.get("last_geneset_update", ""),
                "partial_release_label": analysis.get("partial_release_label", ""),
                "species": analysis.get("species", ""),
                "taxon_id": analysis.get("taxon_id", ""),
                "ref_seqs": analysis.get("ref_seqs", ""),
                "analysis_links": analysis.get("analysis_links", ""),
                "analysis_attributes": analysis.get("analysis_attributes", ""),
                "analysis_type": analysis.get("analysis_type", ""),
                "omit_run_refs_in_test": analysis.get("omit_run_refs_in_test", ""),
                "file_path": file_path,
                "file_type": file_type,
                "remote_name": remote_name,
                "run_accession": run_accession,
                "sample_accession": (row.get("sample_accession") or "").strip(),
                "experiment_accession": (
                    row.get("experiment_accession") or row.get("experiment_accessions") or ""
                ).strip(),
            }
        )

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(args.out)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
