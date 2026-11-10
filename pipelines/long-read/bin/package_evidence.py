#!/usr/bin/env python3
"""Package explicit long-read model/alignment evidence records."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

from parquet_writer import write_string_table


POLICY = {
    "dictionary_only": False,
    "matcher": "explicit_channel_relation",
    "compatibility_policy": "reference_dictionary_exact",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def stable_id(kind: str, payload: dict[str, Any]) -> str:
    encoded = json.dumps({"kind": kind, **payload}, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"invalid JSON evidence input {path}: {error}") from error
    if not isinstance(value, dict):
        raise SystemExit(f"evidence input must be a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def copy_with_checksum(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return digest(destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models-dir", required=True, type=Path)
    parser.add_argument("--alignment-manifests-dir", required=True, type=Path)
    parser.add_argument("--normalized-manifest", required=True, type=Path)
    parser.add_argument("--reference-identity", required=True, type=Path)
    parser.add_argument("--model-records", required=True, type=Path)
    parser.add_argument("--direct-support-dir", type=Path)
    parser.add_argument("--joined-support-dir", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if bool(args.direct_support_dir) != bool(args.joined_support_dir):
        raise SystemExit("direct and joined support must be supplied together")

    records = json.loads(args.model_records.read_text())
    if not isinstance(records, list) or not records:
        raise SystemExit("model-records must be a non-empty JSON array")
    reference_identity = read_json(args.reference_identity)
    reference_sha256 = reference_identity.get("sha256") or digest(args.reference_identity)
    reference_id = stable_id("reference", {"sha256": reference_sha256, "identity": reference_identity})
    evidence = args.output_dir
    evidence.mkdir(parents=True, exist_ok=True)

    reference_path = evidence / "references" / reference_id / "reference.json"
    reference_checksum = copy_with_checksum(args.reference_identity, reference_path)

    normalized_manifest_checksum = digest(args.normalized_manifest)
    shutil.copyfile(args.normalized_manifest, evidence / "normalized_manifest.tsv")

    alignments: dict[str, dict[str, Any]] = {}
    for source in sorted(args.alignment_manifests_dir.glob("*.alignment_manifest.json")):
        payload = read_json(source)
        run_id = str(payload.get("run_id") or source.stem.split(".", 1)[0])
        alignment_id = stable_id("alignment", {"run_id": run_id, "manifest_sha256": digest(source), "reference_id": reference_id})
        dataset_id = run_id
        destination = evidence / "alignments" / alignment_id / dataset_id / "alignment_manifest.json"
        checksum = copy_with_checksum(source, destination)
        alignments[run_id] = {"alignment_id": alignment_id, "dataset_id": dataset_id, "path": str(destination.relative_to(evidence)), "sha256": checksum, "payload": payload}

    staged_models = {path.name: path for path in args.models_dir.glob("*") if path.is_file()}
    models: list[dict[str, Any]] = []
    relations: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            raise SystemExit("each model record must be an object")
        bed12_name = record.get("bed12_name")
        if not bed12_name or bed12_name not in staged_models:
            raise SystemExit(f"model record does not identify a staged BED12 file: {bed12_name}")
        bed12 = staged_models[bed12_name]
        backend = str(record["backend"])
        scope = str(record["scope"])
        run_id = str(record.get("run_id") or record.get("cohort_id") or "unknown")
        dataset_id = run_id if scope == "run" else str(record.get("cohort_id") or run_id)
        model_id = stable_id("model", {"bed12_sha256": digest(bed12), "backend": backend, "scope": scope, "reference_id": reference_id, "dataset_id": dataset_id})
        destination = evidence / "models" / model_id / dataset_id / f"{dataset_id}.{backend}.{scope}.bed12"
        model_checksum = copy_with_checksum(bed12, destination)
        model_entry = {
            "model_id": model_id,
            "dataset_id": dataset_id,
            "backend": backend,
            "scope": scope,
            "run_id": record.get("run_id"),
            "cohort_id": record.get("cohort_id"),
            "path": str(destination.relative_to(evidence)),
            "sha256": model_checksum,
            "model_manifest": record.get("model_manifest"),
            "validation_report": record.get("validation_report"),
        }
        write_json(evidence / "models" / model_id / dataset_id / "model_record.json", model_entry)
        models.append(model_entry)

        matching_alignments = [alignments[record["run_id"]]] if record.get("run_id") in alignments else []
        if scope == "cohort" and not matching_alignments:
            matching_alignments = [alignments[key] for key in sorted(alignments)]
        for alignment in matching_alignments:
            compatibility = {"status": "compatible", "level": "exact", **POLICY}
            relation_payload = {
                "reference_id": reference_id,
                "alignment_id": alignment["alignment_id"],
                "model_id": model_id,
                "dependency_class": "model_alignment",
                "compatibility": compatibility,
            }
            relation_id = stable_id("assessment", relation_payload)
            assessment_dataset_id = stable_id("assessment_dataset", relation_payload)
            assessment = {"relation_id": relation_id, **relation_payload, "dataset_id": assessment_dataset_id, "matcher_policy": POLICY}
            write_json(evidence / "assessments" / relation_id / assessment_dataset_id / "assessment.json", assessment)
            relations.append(assessment)

    relation_fields = ["relation_id", "dataset_id", "reference_id", "alignment_id", "model_id", "dependency_class", "compatibility_status", "compatibility_level", "dictionary_only", "matcher", "compatibility_policy", "matcher_policy"]
    relation_rows = []
    for relation in sorted(relations, key=lambda item: item["relation_id"]):
        compatibility = relation["compatibility"]
        relation_rows.append({
            "relation_id": relation["relation_id"],
            "dataset_id": relation["dataset_id"],
            "reference_id": relation["reference_id"],
            "alignment_id": relation["alignment_id"],
            "model_id": relation["model_id"],
            "dependency_class": relation["dependency_class"],
            "compatibility_status": compatibility["status"],
            "compatibility_level": compatibility["level"],
            "dictionary_only": str(compatibility["dictionary_only"]).lower(),
            "matcher": compatibility["matcher"],
            "compatibility_policy": compatibility["compatibility_policy"],
            "matcher_policy": json.dumps(compatibility, sort_keys=True, separators=(",", ":")),
        })
    with (evidence / "relations.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=relation_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(relation_rows)
    write_string_table(evidence / "relations.parquet", relation_fields, relation_rows)

    catalog = {
        "schema": "long-read-evidence-v1",
        "reference": {"reference_id": reference_id, "path": str(reference_path.relative_to(evidence)), "sha256": reference_checksum},
        "normalized_manifest": {"path": "normalized_manifest.tsv", "sha256": normalized_manifest_checksum},
        "alignments": alignments,
        "models": models,
        "assessments": relation_rows,
        "relations": {"path": "relations.parquet", "inspection_path": "relations.tsv", "rows": len(relation_rows)},
    }
    if args.direct_support_dir:
        support_root = evidence / "support"
        support_root.mkdir()
        files = {}
        run_dirs = sorted(path for path in args.direct_support_dir.iterdir() if path.is_dir())
        if not run_dirs:
            raise SystemExit("direct support contains no run directories")
        for directory in run_dirs:
            metadata = read_json(directory / "evidence_metadata.json")
            if metadata.get("schema") != "direct-alignment-support-v1":
                raise SystemExit(f"unsupported direct support schema: {directory}")
            run_id = metadata["run_id"]
            if run_id not in alignments:
                raise SystemExit(f"direct support has no alignment manifest: {run_id}")
            for source in sorted(directory.iterdir()):
                if not source.is_file():
                    raise SystemExit(f"unexpected direct support entry: {source}")
                relative = Path("per_run") / run_id / source.name
                files[str(relative)] = copy_with_checksum(source, support_root / relative)
        joined_files = sorted(args.joined_support_dir.iterdir())
        required_joined = {"model_index.tsv.gz", "model_support.tsv.gz", "model_junction_support.tsv.gz", "model_event_support.tsv.gz", "chain_event_support.tsv.gz", "event_caller_run_support.tsv.gz", "event_support.tsv.gz", "event_support_sensitivity.tsv.gz", "event_cluster_sensitivity.tsv.gz", "event_tier_summary.tsv", "event_best_depth_summary.tsv", "chain_length_depth_composition.tsv.gz", "support_schema.json"}
        if {path.name for path in joined_files} != required_joined:
            raise SystemExit("joined support files differ from the required contract")
        for source in joined_files:
            relative = Path("summary") / source.name
            files[str(relative)] = copy_with_checksum(source, support_root / relative)
        recording = {
            "schema": "long-read-support-recording-v1",
            "reference_id": reference_id,
            "run_ids": sorted(metadata["run_id"] for metadata in (read_json(path / "evidence_metadata.json") for path in run_dirs)),
            "files": files,
        }
        write_json(support_root / "recording_manifest.json", recording)
        files["recording_manifest.json"] = digest(support_root / "recording_manifest.json")
        (support_root / "SHA256SUMS").write_text("".join(f"{value}  {name}\n" for name, value in sorted(files.items())))
        catalog["direct_support"] = {
            "schema": recording["schema"],
            "path": "support/recording_manifest.json",
            "checksum_path": "support/SHA256SUMS",
            "files": len(files),
        }
    write_json(evidence / "catalog.json", catalog)


if __name__ == "__main__":
    main()
