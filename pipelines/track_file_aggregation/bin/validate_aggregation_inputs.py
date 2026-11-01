#!/usr/bin/env python3
"""Validate aggregation membership and expand it to source-track rows."""

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
    "tool_versions",
    "normalization_parameters",
}
AGGREGATION_FIELDS = {
    "aggregate_id",
    "aggregate_type",
    "source_entity_id",
    "source_entity_type",
    "aggregation_method",
}
AGGREGATE_TYPES = {"sample", "tissue", "global"}
AGGREGATION_METHODS = {"raw", "equal_sample"}
TRACK_TYPES = {
    "coverage_forward",
    "coverage_reverse",
    "splice_junction",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        return list(reader)


def canonical_json(raw: str, field: str) -> str:
    if not raw.strip():
        return "{}"
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{field} is not valid JSON: {raw!r}") from exc
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def normalization_signature(raw: str, track_type: str, allow_opaque: bool) -> str:
    if not track_type.startswith("coverage_"):
        return "not_applicable"
    if not raw.strip():
        raise ValueError("coverage track has empty normalization_parameters")
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        if not allow_opaque:
            raise ValueError(
                "coverage normalization_parameters must be JSON; "
                "regenerate the track manifest or set allow_opaque_normalization=true "
                "for a controlled compatibility test"
            )
        return f"opaque:{raw}"
    if not isinstance(value, dict):
        raise ValueError("coverage normalization_parameters must contain a JSON object")

    # Aggregate manifests retain the source preparation metadata under
    # source_normalization_parameters so sample-level tracks can be consumed by
    # a later tissue-level aggregation.
    source_value = value
    if "bin_size" not in source_value or "normalization" not in source_value:
        source_value = value.get("source_normalization_parameters")
        if isinstance(source_value, str):
            try:
                source_value = json.loads(source_value)
            except json.JSONDecodeError:
                if not allow_opaque:
                    raise ValueError(
                        "coverage source normalization_parameters must be JSON; "
                        "regenerate the track manifest or set allow_opaque_normalization=true "
                        "for a controlled compatibility test"
                    )
        if isinstance(source_value, str) and source_value.startswith("opaque:"):
            return source_value
    if not isinstance(source_value, dict):
        if allow_opaque:
            return f"opaque:{source_value}"
        raise ValueError("coverage normalization_parameters must contain a JSON object")
    comparable = {key: item for key, item in source_value.items() if key not in {"strand"}}
    if "bin_size" not in comparable or "normalization" not in comparable:
        if allow_opaque and isinstance(value.get("source_normalization_parameters"), str):
            return f"opaque:{value['source_normalization_parameters']}"
        raise ValueError("coverage normalization_parameters must include bin_size and normalization")
    return json.dumps(comparable, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def safe_id(value: str) -> str:
    result = "".join(char if char.isalnum() or char in "_.-" else "_" for char in value)
    return result or "aggregate"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def build_plan(
    track_rows: list[dict[str, str]],
    aggregation_rows: list[dict[str, str]],
    base_dir: Path,
    selected: set[str],
    allow_opaque: bool,
) -> tuple[list[dict[str, str]], list[dict[str, object]], list[str]]:
    missing = TRACK_FIELDS - set(track_rows[0]) if track_rows else TRACK_FIELDS
    if missing:
        raise ValueError(f"track manifest is missing columns: {', '.join(sorted(missing))}")
    missing = AGGREGATION_FIELDS - set(aggregation_rows[0]) if aggregation_rows else AGGREGATION_FIELDS
    if missing:
        raise ValueError(f"aggregation manifest is missing columns: {', '.join(sorted(missing))}")

    tracks: dict[tuple[str, str, str, str], dict[str, str]] = {}
    assembly_releases: defaultdict[str, set[str]] = defaultdict(set)
    for row in track_rows:
        key = (row["gca_accession"], row["entity_id"], row["entity_type"], row["track_type"])
        if key in tracks:
            raise ValueError(f"duplicate track manifest row: {key}")
        if row["track_type"] in TRACK_TYPES and row["status"] != "complete":
            raise ValueError(f"source track is not complete: {key} ({row['status']})")
        if row["track_type"] in TRACK_TYPES:
            if not row["sha256"].strip():
                raise ValueError(f"source track has no checksum: {key}")
            if not row["assembly_release"].strip():
                raise ValueError(f"source track has no assembly release: {key}")
            assembly_releases[row["gca_accession"]].add(row["assembly_release"])
            track_path = Path(row["track_path"])
            if not track_path.is_absolute():
                track_path = (base_dir / track_path).resolve()
            if not track_path.is_file():
                raise ValueError(f"source track does not exist: {track_path}")
            row = dict(row)
            row["track_path"] = str(track_path)
        tracks[key] = row
    for gca, releases in assembly_releases.items():
        if len(releases) != 1:
            raise ValueError(f"track manifest has multiple assembly releases for {gca}: {sorted(releases)}")

    group_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    group_config: dict[tuple[str, str], dict[str, str]] = {}
    seen_memberships: set[tuple[str, str, str]] = set()
    assemblies: set[str] = set()
    for row in aggregation_rows:
        aggregate_id = row["aggregate_id"].strip()
        aggregate_type = row["aggregate_type"].strip()
        source_id = row["source_entity_id"].strip()
        source_type = row["source_entity_type"].strip()
        aggregation_method = row["aggregation_method"].strip()
        if (
            not aggregate_id
            or aggregate_type not in AGGREGATE_TYPES
            or not source_id
            or not source_type
            or aggregation_method not in AGGREGATION_METHODS
        ):
            raise ValueError(f"invalid aggregation membership row: {row}")
        if aggregation_method == "equal_sample" and source_type != "sample":
            raise ValueError(
                "equal_sample aggregation requires source_entity_type=sample; "
                "aggregate runs into samples with raw first"
            )
        membership_key = (aggregate_id, source_id, source_type)
        if membership_key in seen_memberships:
            raise ValueError(f"duplicate aggregation membership: {membership_key}")
        seen_memberships.add(membership_key)

        candidate_keys = [key for key in tracks if key[1] == source_id and key[2] == source_type]
        if not candidate_keys:
            raise ValueError(f"source entity is absent from track manifest: {source_type}:{source_id}")
        gca_values = {key[0] for key in candidate_keys}
        if len(gca_values) != 1:
            raise ValueError(f"source entity resolves to multiple assemblies: {source_type}:{source_id}")
        gca = next(iter(gca_values))
        declared_gca = row.get("gca_accession", "").strip()
        if declared_gca and declared_gca != gca:
            raise ValueError(f"assembly mismatch for {source_type}:{source_id}: {declared_gca} != {gca}")
        assemblies.add(gca)

        group_key = (gca, aggregate_id)
        config = {
            "aggregate_type": aggregate_type,
            "aggregate_label": row.get("aggregate_label", "").strip() or aggregate_id,
            "aggregation_method": aggregation_method,
            "metadata_json": canonical_json(row.get("metadata_json", ""), "metadata_json"),
        }
        if group_key in group_config and group_config[group_key] != config:
            raise ValueError(f"aggregate metadata or policy differs within group: {group_key}")
        group_config[group_key] = config
        group_rows[group_key].append({"source_id": source_id, "source_type": source_type, "gca": gca})

    if len(assemblies) != 1:
        raise ValueError("one assembly is required per aggregation invocation")

    aggregate_safe_ids: dict[tuple[str, str], str] = {}
    for gca, aggregate_id in group_rows:
        safe_aggregate_id = safe_id(aggregate_id)
        for (existing_gca, existing_id), existing_safe_id in aggregate_safe_ids.items():
            if existing_gca == gca and existing_safe_id == safe_aggregate_id and existing_id != aggregate_id:
                raise ValueError(
                    f"aggregate IDs collide after sanitization for {gca}: "
                    f"{existing_id!r} and {aggregate_id!r} -> {safe_aggregate_id!r}"
                )
        aggregate_safe_ids[(gca, aggregate_id)] = safe_aggregate_id

    plan: list[dict[str, str]] = []
    metadata: list[dict[str, object]] = []
    for group_key in sorted(group_rows):
        gca, aggregate_id = group_key
        members = group_rows[group_key]
        config = group_config[group_key]
        source_ids = [member["source_id"] for member in members]
        source_types = {member["source_type"] for member in members}
        if len(source_types) != 1:
            raise ValueError(f"aggregate mixes source entity types: {group_key}")
        entity_type = next(iter(source_types))
        selected_rows: list[dict[str, str]] = []
        for member in members:
            source_tracks = {
                track_type: tracks[(gca, member["source_id"], entity_type, track_type)]
                for track_type in TRACK_TYPES
                if (gca, member["source_id"], entity_type, track_type) in tracks
            }
            if "coverage" in selected and not {"coverage_forward", "coverage_reverse"}.issubset(source_tracks):
                raise ValueError(f"source lacks a complete forward/reverse coverage pair: {member['source_id']}")
            if "splice_junction" in selected and "splice_junction" not in source_tracks:
                raise ValueError(f"source lacks a splice-junction track: {member['source_id']}")
            selected_rows.extend(
                track
                for track_type, track in source_tracks.items()
                if (track_type.startswith("coverage_") and "coverage" in selected)
                or (track_type == "splice_junction" and "splice_junction" in selected)
            )

        coverage_signatures = {
            normalization_signature(row["normalization_parameters"], row["track_type"], allow_opaque)
            for row in selected_rows
            if row["track_type"].startswith("coverage_") and "coverage" in selected
        }
        if len(coverage_signatures) > 1:
            raise ValueError(f"coverage normalization differs within aggregate: {group_key}")
        for row in selected_rows:
            if row["track_type"].startswith("coverage_"):
                signature = normalization_signature(row["normalization_parameters"], row["track_type"], allow_opaque)
            else:
                signature = "not_applicable"
            plan.append(
                {
                    "gca_accession": gca,
                    "assembly_release": row["assembly_release"],
                    "aggregate_id": aggregate_id,
                    "aggregate_type": config["aggregate_type"],
                    "aggregate_label": config["aggregate_label"],
                    "source_entity_id": row["entity_id"],
                    "source_entity_type": row["entity_type"],
                    "track_type": row["track_type"],
                    "track_path": row["track_path"],
                    "source_path": row["source_path"],
                    "expected_sha256": row["sha256"],
                    "normalization_parameters": quote(row["normalization_parameters"], safe=""),
                    "normalization_signature": quote(signature, safe=""),
                    "aggregation_method": config["aggregation_method"],
                    "metadata_json": quote(config["metadata_json"], safe=""),
                }
            )
        metadata.append(
            {
                "gca_accession": gca,
                "aggregate_id": aggregate_id,
                "aggregate_type": config["aggregate_type"],
                "aggregate_label": config["aggregate_label"],
                "source_entity_ids": source_ids,
                "source_entity_type": entity_type,
                "source_count": len(source_ids),
                "aggregation_method": config["aggregation_method"],
                "metadata": json.loads(config["metadata_json"]),
            }
        )
    if not plan:
        raise ValueError("aggregation manifest has no source tracks for the selected track_types")
    return plan, metadata, sorted(assemblies)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--track-manifest", type=Path, required=True)
    parser.add_argument("--aggregation-manifest", type=Path, required=True)
    parser.add_argument("--manifest-base-dir", type=Path, required=True)
    parser.add_argument("--track-types", choices=["all", "coverage", "splice_junction"], default="all")
    parser.add_argument("--allow-opaque-normalization", action="store_true")
    parser.add_argument("--output-plan", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--output-validation", type=Path, required=True)
    args = parser.parse_args()

    selected = {"coverage", "splice_junction"} if args.track_types == "all" else {args.track_types}
    plan, metadata, assemblies = build_plan(
        read_tsv(args.track_manifest),
        read_tsv(args.aggregation_manifest),
        args.manifest_base_dir.resolve(),
        selected,
        args.allow_opaque_normalization,
    )
    fields = list(plan[0]) if plan else [
        "gca_accession", "assembly_release", "aggregate_id", "aggregate_type", "aggregate_label",
        "source_entity_id", "source_entity_type", "track_type", "track_path", "expected_sha256",
        "source_path", "normalization_parameters", "normalization_signature", "aggregation_method", "metadata_json",
    ]
    write_tsv(args.output_plan, fields, plan)
    args.output_metadata.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_tsv(
        args.output_validation,
        ["status", "assembly", "source_tracks", "aggregate_groups"],
        [{"status": "complete", "assembly": assemblies[0], "source_tracks": str(len(plan)), "aggregate_groups": str(len(metadata))}],
    )


if __name__ == "__main__":
    main()
