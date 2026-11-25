#!/usr/bin/env python3
"""Validate a run-level track manifest and write a QC aggregation plan."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote


TRACK_FIELDS = {
    "gca_accession",
    "assembly_release",
    "entity_id",
    "entity_type",
    "track_type",
    "track_path",
    "source_path",
    "status",
    "sha256",
    "normalization_parameters",
}
TRACK_TYPES = {"coverage_forward", "coverage_reverse", "splice_junction"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path} has no header")
        return list(reader)


def safe_id(value: str) -> str:
    result = "".join(char if char.isalnum() or char in "_.-" else "_" for char in value)
    return result or "run"


def normalization_signature(raw: str, allow_opaque: bool) -> str:
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        if allow_opaque:
            return f"opaque:{raw}"
        raise ValueError("coverage normalization_parameters must be valid JSON")
    if not isinstance(value, dict):
        raise ValueError("coverage normalization_parameters must be a JSON object")
    if "bin_size" not in value or "normalization" not in value:
        if allow_opaque:
            return f"opaque:{raw}"
        raise ValueError("coverage normalization_parameters must include bin_size and normalization")
    comparable = {key: item for key, item in value.items() if key != "strand"}
    return json.dumps(comparable, sort_keys=True, separators=(",", ":"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def build_plan(
    rows: list[dict[str, str]],
    base_dir: Path,
    selected: set[str],
    entity_type: str,
    species_name: str,
    allow_opaque: bool,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    if not rows:
        raise ValueError("track manifest has no rows")
    missing = TRACK_FIELDS - set(rows[0])
    if missing:
        raise ValueError(f"track manifest is missing columns: {', '.join(sorted(missing))}")

    selected_track_types = {
        track_type
        for track_type in TRACK_TYPES
        if (track_type.startswith("coverage_") and "coverage" in selected)
        or (track_type == "splice_junction" and "splice_junction" in selected)
    }
    tracks: dict[tuple[str, str], dict[str, str]] = {}
    assemblies: set[tuple[str, str]] = set()
    safe_ids: dict[str, str] = {}
    for row in rows:
        if row["track_type"] not in selected_track_types:
            continue
        key = (row["entity_id"], row["track_type"])
        if key in tracks:
            raise ValueError(f"duplicate track manifest row: {key}")
        if row["entity_type"] != entity_type:
            continue
        if row["status"] != "complete":
            raise ValueError(f"source track is not complete: {key} ({row['status']})")
        if not row["sha256"].strip():
            raise ValueError(f"source track has no checksum: {key}")
        path = Path(row["track_path"])
        if not path.is_absolute():
            path = (base_dir / path).resolve()
        if not path.is_file():
            raise ValueError(f"source track does not exist: {path}")
        safe = safe_id(row["entity_id"])
        previous = safe_ids.get(safe)
        if previous and previous != row["entity_id"]:
            raise ValueError(f"entity IDs collide after sanitization: {previous!r} and {row['entity_id']!r}")
        safe_ids[safe] = row["entity_id"]
        if row["track_type"].startswith("coverage_"):
            normalization_signature(row["normalization_parameters"], allow_opaque)
        assemblies.add((row["gca_accession"], row["assembly_release"]))
        item = dict(row)
        item["track_path"] = str(path)
        item["safe_id"] = safe
        tracks[key] = item

    if len(assemblies) != 1:
        raise ValueError("one assembly and assembly release are required per invocation")
    if not tracks:
        raise ValueError("no complete source tracks match the selected entity and track types")

    source_ids = sorted({entity_id for entity_id, _track_type in tracks})
    for source_id in source_ids:
        if "coverage" in selected and (source_id, "coverage_forward") not in tracks:
            raise ValueError(f"source entity lacks coverage_forward: {source_id}")
        if "coverage" in selected and (source_id, "coverage_reverse") not in tracks:
            raise ValueError(f"source entity lacks coverage_reverse: {source_id}")
        if "splice_junction" in selected and (source_id, "splice_junction") not in tracks:
            raise ValueError(f"source entity lacks splice_junction: {source_id}")

    coverage_signatures = {
        normalization_signature(item["normalization_parameters"], allow_opaque)
        for item in tracks.values()
        if item["track_type"].startswith("coverage_")
    }
    if len(coverage_signatures) > 1:
        raise ValueError("coverage normalization differs between source tracks")

    plan = []
    for key in sorted(tracks):
        item = tracks[key]
        plan.append({
            "gca_accession": item["gca_accession"],
            "assembly_release": item["assembly_release"],
            "species_name": species_name or item["gca_accession"],
            "entity_id": item["entity_id"],
            "entity_type": item["entity_type"],
            "safe_id": item["safe_id"],
            "track_type": item["track_type"],
            "track_path": item["track_path"],
            "source_path": item["source_path"],
            "sha256": item["sha256"],
            "normalization_parameters": quote(item["normalization_parameters"], safe=""),
        })
    gca, release = next(iter(assemblies))
    metadata = {
        "gca_accession": gca,
        "assembly_release": release,
        "species_name": species_name or gca,
        "entity_type": entity_type,
        "source_entity_ids": source_ids,
        "source_count": len(source_ids),
        "track_types": sorted(selected),
    }
    return plan, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--track-manifest", type=Path, required=True)
    parser.add_argument("--manifest-base-dir", type=Path, required=True)
    parser.add_argument("--track-types", choices=["all", "coverage", "splice_junction"], default="all")
    parser.add_argument("--entity-type", default="run")
    parser.add_argument("--species-name", default="")
    parser.add_argument("--allow-opaque-normalization", action="store_true")
    parser.add_argument("--output-plan", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--output-validation", type=Path, required=True)
    args = parser.parse_args()

    selected = {"coverage", "splice_junction"} if args.track_types == "all" else {args.track_types}
    plan, metadata = build_plan(
        read_tsv(args.track_manifest),
        args.manifest_base_dir.resolve(),
        selected,
        args.entity_type,
        args.species_name,
        args.allow_opaque_normalization,
    )
    fields = list(plan[0])
    write_tsv(args.output_plan, fields, plan)
    args.output_metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    write_tsv(
        args.output_validation,
        ["status", "gca_accession", "assembly_release", "source_tracks", "source_entities"],
        [{
            "status": "complete",
            "gca_accession": metadata["gca_accession"],
            "assembly_release": metadata["assembly_release"],
            "source_tracks": str(len(plan)),
            "source_entities": str(metadata["source_count"]),
        }],
    )


if __name__ == "__main__":
    main()
