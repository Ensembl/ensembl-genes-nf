#!/usr/bin/env python3
"""Write output inventory and provenance for the custom QC evidence bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_plan(path: Path) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("supporting evidence plan is empty")
    return {
        "gca_accession": rows[0]["gca_accession"],
        "assembly_release": rows[0]["assembly_release"],
        "species_name": rows[0].get("species_name") or rows[0]["gca_accession"],
        "source_entity_ids": sorted({row["entity_id"] for row in rows}),
        "source_track_count": len(rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--files", nargs="+", type=Path, required=True)
    args = parser.parse_args()

    plan = read_plan(args.plan)
    output_files = [path for path in args.files if not path.name.endswith(".descriptor.tsv")]
    rows = []
    for path in output_files:
        kind = "coverage" if path.name.endswith(".npz") or ".coverage_" in path.name else "splice_junction"
        rows.append({
            "gca_accession": plan["gca_accession"],
            "assembly_release": plan["assembly_release"],
            "species_name": plan["species_name"],
            "output_kind": kind,
            "output_path": path.name,
            "source_entity_count": str(len(plan["source_entity_ids"])),
            "source_entity_ids": ",".join(plan["source_entity_ids"]),
            "status": "complete",
            "sha256": sha256(path),
        })

    fields = list(rows[0]) if rows else [
            "gca_accession", "assembly_release", "species_name", "output_kind", "output_path",
        "source_entity_count", "source_entity_ids", "status", "sha256",
    ]
    with Path("supporting_evidence_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with Path("supporting_evidence_status.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gca_accession", "assembly_release", "species_name", "status", "source_entity_count", "output_count"])
        writer.writerow([
            plan["gca_accession"], plan["assembly_release"], plan["species_name"], "complete",
            len(plan["source_entity_ids"]), len(rows),
        ])
    Path("supporting_evidence_metadata.json").write_text(json.dumps({
        **plan,
        "outputs": rows,
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
