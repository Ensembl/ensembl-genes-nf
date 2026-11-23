#!/usr/bin/env python3
"""Create deterministic, diversity-first production manifest tiers."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


UNKNOWN = "<unknown>"
REQUIRED_COLUMNS = {"run_id", "technology"}


def normalise(value: str | None) -> str:
    value = (value or "").strip().casefold()
    value = re.sub(r"\s+", " ", value)
    return value or UNKNOWN


def biological_context(row: dict[str, str]) -> str:
    for field in ("tissue", "organism_part", "cell_type", "developmental_stage"):
        value = normalise(row.get(field))
        if value != UNKNOWN:
            return f"{field}:{value}"
    return UNKNOWN


def features(row: dict[str, str]) -> dict[str, str]:
    return {
        "tissue": biological_context(row),
        "technology": normalise(row.get("technology")),
        "study": normalise(row.get("study_accession")),
        "sample": normalise(row.get("sample_id")),
        "strain_or_breed": normalise(row.get("strain_or_breed")),
        "read_representation": normalise(row.get("read_representation")),
    }


def score(
    row: dict[str, str],
    selected: list[dict[str, str]],
    technology_targets: dict[str, int],
    quota_selected: list[dict[str, str]],
) -> tuple[int, int]:
    current = features(row)
    seen = {field: {features(item)[field] for item in selected} for field in current}
    weights = {
        "tissue": 100,
        "technology": 30,
        "study": 20,
        "sample": 10,
        "strain_or_breed": 5,
        "read_representation": 5,
    }
    novelty = sum(
        weight
        for field, weight in weights.items()
        if current[field] != UNKNOWN and current[field] not in seen[field]
    )
    technology_count = sum(
        1 for item in quota_selected if features(item)["technology"] == current["technology"]
    )
    quota_bonus = 250 if technology_count < technology_targets.get(current["technology"], 0) else 0
    known_biology = int(current["tissue"] != UNKNOWN)
    return novelty + quota_bonus, known_biology


def choose(
    rows: list[dict[str, str]],
    count: int,
    context: list[dict[str, str]] | None = None,
    technology_targets: dict[str, int] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, int | str]]]:
    selected: list[dict[str, str]] = []
    audit: list[dict[str, int | str]] = []
    remaining = list(rows)
    selected_context = list(context or [])
    quota_selected: list[dict[str, str]] = []
    technology_targets = technology_targets or {}

    while remaining and len(selected) < count:
        ranked = sorted(
            remaining,
            key=lambda row: (
                -score(row, selected_context, technology_targets, quota_selected)[0],
                -score(row, selected_context, technology_targets, quota_selected)[1],
                normalise(row.get("run_id")),
            ),
        )
        row = ranked[0]
        novelty, known_biology = score(row, selected_context, technology_targets, quota_selected)
        feature_context = features(row)
        selected.append(row)
        selected_context.append(row)
        quota_selected.append(row)
        remaining.remove(row)
        audit.append(
            {
                "run_id": row["run_id"],
                "selection_rank": len(selected),
                "selection_score": novelty,
                "known_biology": known_biology,
                "tissue_context": feature_context["tissue"],
                "technology": feature_context["technology"],
                "study_accession": feature_context["study"],
                "sample_id": feature_context["sample"],
            }
        )
    return selected, audit


def read_targets(path: Path) -> dict[str, tuple[int, int]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"species", "tier1_count", "tier2_count", "tier1_pacbio", "tier2_pacbio"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} must contain species, tier1_count, and tier2_count columns")
    targets = {}
    for row in rows:
        species = row["species"].strip()
        tier1 = int(row["tier1_count"])
        tier2 = int(row["tier2_count"])
        tier1_pacbio = int(row["tier1_pacbio"])
        tier2_pacbio = int(row["tier2_pacbio"])
        if min(tier1, tier2, tier1_pacbio, tier2_pacbio) < 0:
            raise ValueError(f"Negative tier target for {species}")
        if tier1_pacbio > tier1 or tier2_pacbio > tier2:
            raise ValueError(f"PacBio target exceeds tier size for {species}")
        targets[species] = (tier1, tier2, tier1_pacbio, tier2_pacbio)
    return targets


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-root", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    args = parser.parse_args()

    targets = read_targets(args.targets)
    index_rows = []
    audit_rows = []

    for species in sorted(targets):
        primary = args.manifest_root / species / "primary.tsv"
        if not primary.exists():
            raise FileNotFoundError(primary)
        with primary.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            missing = REQUIRED_COLUMNS.difference(fieldnames)
            if missing:
                raise ValueError(f"{primary} is missing columns: {sorted(missing)}")
            rows = list(reader)

        tier1_target, tier2_target, tier1_pacbio, tier2_pacbio = targets[species]
        tier1, tier1_audit = choose(
            rows,
            min(tier1_target, len(rows)),
            technology_targets={"pacbio": tier1_pacbio, "ont": tier1_target - tier1_pacbio},
        )
        tier1_ids = {row["run_id"] for row in tier1}
        after_tier1 = [row for row in rows if row["run_id"] not in tier1_ids]
        tier2, tier2_audit = choose(
            after_tier1,
            min(tier2_target, len(after_tier1)),
            context=tier1,
            technology_targets={"pacbio": tier2_pacbio, "ont": tier2_target - tier2_pacbio},
        )
        selected_ids = {row["run_id"] for row in tier1 + tier2}
        remaining = [row for row in rows if row["run_id"] not in selected_ids]

        write_rows(primary.parent / "tier1.tsv", fieldnames, tier1)
        write_rows(primary.parent / "tier2.tsv", fieldnames, tier2)
        write_rows(primary.parent / "remaining.tsv", fieldnames, remaining)

        for tier, audit in ((1, tier1_audit), (2, tier2_audit)):
            for item in audit:
                audit_rows.append({"species": species, "tier": tier, **item})

        index_rows.append(
            {
                "species": species,
                "primary_runs": len(rows),
                "tier1_runs": len(tier1),
                "tier2_runs": len(tier2),
                "remaining_runs": len(remaining),
                "tier1_manifest": f"{species}/tier1.tsv",
                "tier2_manifest": f"{species}/tier2.tsv",
                "remaining_manifest": f"{species}/remaining.tsv",
            }
        )

    write_rows(
        args.manifest_root / "tier_manifest_index.tsv",
        list(index_rows[0]),
        index_rows,
    )
    write_rows(
        args.manifest_root / "tier_selection_audit.tsv",
        [
            "species",
            "tier",
            "selection_rank",
            "run_id",
            "selection_score",
            "known_biology",
            "tissue_context",
            "technology",
            "study_accession",
            "sample_id",
        ],
        audit_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
