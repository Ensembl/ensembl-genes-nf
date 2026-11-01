#!/usr/bin/env python3
"""Prepare explicit staging metadata for retrospective evidence collection."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--models-dir", required=True, type=Path)
    parser.add_argument("--alignment-manifests-dir", required=True, type=Path)
    parser.add_argument("--reference-identity", required=True, type=Path)
    parser.add_argument("--model-records", required=True, type=Path)
    args = parser.parse_args()

    if not args.reference.is_file():
        raise SystemExit(f"retrospective reference does not exist: {args.reference}")
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("retrospective manifest contains no rows")

    model_columns = ("bed12", "model_bed12", "model", "bed12_path")
    args.models_dir.mkdir(parents=True, exist_ok=True)
    args.alignment_manifests_dir.mkdir(parents=True, exist_ok=True)
    reference_sha256 = digest(args.reference)
    reference_identity = {
        "path": str(args.reference),
        "sha256": reference_sha256,
        "dictionary_sha256": None,
        "retrospective": True,
    }
    args.reference_identity.write_text(json.dumps(reference_identity, indent=2, sort_keys=True) + "\n")

    records = []
    for row in rows:
        run_id = row.get("run_id") or row.get("id")
        if not run_id:
            raise SystemExit("retrospective rows require run_id or id")
        alignment = row.get("alignment") or row.get("bam") or row.get("input")
        if not alignment or not Path(alignment).is_file():
            raise SystemExit(f"retrospective row {run_id} requires an existing BAM/alignment input")
        alignment_path = Path(alignment)
        alignment_manifest = args.alignment_manifests_dir / f"{run_id}.alignment_manifest.json"
        alignment_manifest.write_text(json.dumps({
            "run_id": run_id,
            "cohort_id": row.get("cohort_id") or args.cohort_id,
            "bam": str(alignment_path),
            "bam_sha256": digest(alignment_path),
            "reference_id": reference_sha256,
            "retrospective": True,
        }, indent=2, sort_keys=True) + "\n")

        model_column = next((name for name in model_columns if row.get(name)), None)
        if model_column is None:
            raise SystemExit(f"retrospective row {run_id} requires one of: {', '.join(model_columns)}")
        model_path = Path(row[model_column])
        if not model_path.is_file():
            raise SystemExit(f"retrospective model artifact does not exist: {model_path}")
        backend = row.get("backend") or "retrospective"
        scope = row.get("scope") or "run"
        staged_name = f"{run_id}.{backend}.{scope}.bed12"
        shutil.copyfile(model_path, args.models_dir / staged_name)
        records.append({
            "bed12_name": staged_name,
            "backend": backend,
            "scope": scope,
            "run_id": run_id,
            "cohort_id": row.get("cohort_id") or args.cohort_id,
            "model_manifest": json.dumps({"source": str(model_path), "retrospective": True}, sort_keys=True),
            "validation_report": "retrospective artifact; validation delegated to source manifest",
        })

    args.model_records.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
