#!/usr/bin/env python3
"""Write one aggregate track record and provenance sidecar."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = [
    "gca_accession", "assembly_release", "entity_id", "entity_type", "track_type", "track_path",
    "source_path", "status", "sha256", "tool_versions", "normalization_parameters", "aggregate_type",
    "aggregate_label", "source_entity_ids", "source_count", "aggregation_method", "metadata_json",
    "provenance_path",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--track", type=Path, required=True)
    parser.add_argument("--checksum", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    meta = json.loads(args.metadata.read_text(encoding="utf-8"))
    checksum = args.checksum.read_text(encoding="utf-8").split()[0]
    track_path = meta["published_track_path"]
    record = {
        "gca_accession": meta["gca_accession"],
        "assembly_release": meta["assembly_release"],
        "entity_id": meta["aggregate_id"],
        "entity_type": meta["aggregate_type"],
        "track_type": meta["current_track"],
        "track_path": track_path,
        "source_path": meta["source_path"],
        "status": "complete",
        "sha256": checksum,
        "tool_versions": "aggregate_pipeline",
        "normalization_parameters": meta["normalization_parameters"],
        "aggregate_type": meta["aggregate_type"],
        "aggregate_label": meta["aggregate_label"],
        "source_entity_ids": meta["source_entity_ids"],
        "source_count": str(meta["source_count"]),
        "aggregation_method": meta["aggregation_method"],
        "metadata_json": meta["metadata_json"],
        "provenance_path": meta["provenance_path"],
    }
    with args.result.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(record)
    args.provenance.write_text(json.dumps({"record": record, "track_sha256": sha256(args.track)}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
