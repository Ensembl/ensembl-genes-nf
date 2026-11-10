#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", required=True, type=Path)
    parser.add_argument("--require-support", action="store_true")
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
        if digest(model_path) != model.get("sha256"):
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
    with (args.evidence_dir / "relations.parquet").open("rb") as handle:
        first = handle.read(4)
        handle.seek(-4, 2)
        last = handle.read(4)
    if first != b"PAR1" or last != b"PAR1":
        raise SystemExit("relations.parquet is not a Parquet file")
    support = catalog.get("direct_support")
    if args.require_support and not support:
        raise SystemExit("direct support is required but absent")
    if support:
        root = args.evidence_dir / "support"
        recording = json.loads((root / "recording_manifest.json").read_text())
        if recording.get("schema") != "long-read-support-recording-v1" or recording.get("reference_id") != reference.get("reference_id"):
            raise SystemExit("direct support recording identity or schema mismatch")
        expected = dict(recording["files"])
        expected["recording_manifest.json"] = digest(root / "recording_manifest.json")
        required_summary = {
            "model_index.tsv.gz", "model_support.tsv.gz", "model_junction_support.tsv.gz",
            "model_event_support.tsv.gz", "chain_event_support.tsv.gz", "event_caller_run_support.tsv.gz",
            "event_support.tsv.gz", "event_support_sensitivity.tsv.gz",
            "event_cluster_sensitivity.tsv.gz", "event_tier_summary.tsv",
            "event_best_depth_summary.tsv", "chain_length_depth_composition.tsv.gz",
            "support_schema.json",
        }
        required_per_run = {
            "alignment_filter_flow.tsv", "alignment_noise.tsv", "chain_observations.tsv.gz",
            "junction_observations.tsv.gz", "end_observations.tsv.gz", "evidence_metadata.json",
        }
        required_paths = {f"summary/{name}" for name in required_summary}
        required_paths.update(f"per_run/{run_id}/{name}" for run_id in recording["run_ids"] for name in required_per_run)
        if set(recording["files"]) != required_paths or support.get("files") != len(expected):
            raise SystemExit("direct support file inventory does not match the contract")
        checksum_rows = (root / "SHA256SUMS").read_text().splitlines()
        if checksum_rows != [f"{value}  {name}" for name, value in sorted(expected.items())]:
            raise SystemExit("direct support checksum list mismatch")
        for name, checksum in expected.items():
            path = root / name
            if not path.is_file() or digest(path) != checksum:
                raise SystemExit(f"direct support checksum mismatch: {name}")
        if set(recording["run_ids"]) != set(alignments):
            raise SystemExit("direct support run coverage does not match alignments")
    payload = {name: digest(args.evidence_dir / name) for name in required}
    verification = {"valid": True, "schema": catalog["schema"], "models": len(models), "relations": len(rows), "checksums": payload}
    (args.evidence_dir / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")


if __name__ == "__main__":
    main()
