#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", required=True, type=Path)
    args = parser.parse_args()
    required = ["catalog.json", "relations.parquet", "relations.tsv"]
    missing = [name for name in required if not (args.evidence_dir / name).is_file()]
    if missing:
        raise SystemExit("missing evidence products: " + ", ".join(missing))
    try:
        catalog = json.loads((args.evidence_dir / "catalog.json").read_text())
    except json.JSONDecodeError as error:
        raise SystemExit(f"invalid catalog.json: {error}") from error
    if catalog.get("schema") != "long-read-evidence-v1":
        raise SystemExit("catalog.json has an unsupported evidence schema")
    reference = catalog.get("reference", {})
    if not (args.evidence_dir / reference.get("path", "")).is_file():
        raise SystemExit("catalog reference artifact is missing")
    models = catalog.get("models", [])
    model_ids = {model.get("model_id") for model in models}
    if len(model_ids) != len(models) or None in model_ids:
        raise SystemExit("catalog model IDs are not unique")
    for model in models:
        model_path = args.evidence_dir / model.get("path", "")
        if not model_path.is_file():
            raise SystemExit(f"catalog model artifact is missing: {model_path}")
        if hashlib.sha256(model_path.read_bytes()).hexdigest() != model.get("sha256"):
            raise SystemExit(f"catalog model checksum mismatch: {model_path}")
    alignments = catalog.get("alignments", {})
    alignment_ids = {value.get("alignment_id") for value in alignments.values()}
    if len(alignment_ids) != len(alignments) or None in alignment_ids:
        raise SystemExit("catalog alignment IDs are not unique")
    for value in alignments.values():
        if not (args.evidence_dir / value.get("path", "")).is_file():
            raise SystemExit("catalog alignment artifact is missing")
    with (args.evidence_dir / "relations.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    relation_ids = {row.get("relation_id") for row in rows}
    if len(relation_ids) != len(rows) or None in relation_ids:
        raise SystemExit("evidence relation IDs are not unique")
    for row in rows:
        if row.get("model_id") not in model_ids or row.get("alignment_id") not in alignment_ids:
            raise SystemExit("relation references an unknown model or alignment ID")
        if not row.get("dataset_id") or len(row["dataset_id"]) != 64 or not row.get("dependency_class") or not row.get("compatibility_level") or not row.get("matcher_policy"):
            raise SystemExit("relation is missing dependency or compatibility metadata")
    parquet_bytes = (args.evidence_dir / "relations.parquet").read_bytes()
    if not parquet_bytes.startswith(b"PAR1") or not parquet_bytes.endswith(b"PAR1"):
        raise SystemExit("relations.parquet is not a Parquet file")
    payload = {name: hashlib.sha256((args.evidence_dir / name).read_bytes()).hexdigest() for name in required}
    verification = {"valid": True, "schema": catalog["schema"], "models": len(models), "relations": len(rows), "checksums": payload}
    (args.evidence_dir / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")


if __name__ == "__main__":
    main()
