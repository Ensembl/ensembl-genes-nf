#!/usr/bin/env python3
"""Write the final run manifest and staged model checksum table."""
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
    parser.add_argument("--configuration", required=True)
    parser.add_argument("--stage-status", required=True)
    parser.add_argument("--model-publications", required=True)
    parser.add_argument("--evidence-status", required=True)
    parser.add_argument("--merge-provenance", required=True)
    parser.add_argument("--normalized-manifest", required=True, type=Path)
    parser.add_argument("--input-checksums", required=True, type=Path)
    parser.add_argument("--preparation-status", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--reference-identity", required=True, type=Path)
    parser.add_argument("--model-reports", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--checksums", required=True, type=Path)
    args = parser.parse_args()

    status = json.loads(args.stage_status)
    status["evidence"] = args.evidence_status
    status["finalization"] = "COMPLETE"
    config = json.loads(args.configuration)
    model_publications = json.loads(args.model_publications)
    merge_provenance = json.loads(args.merge_provenance)
    with args.normalized_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    run_ids = sorted({row.get("run_id") for row in rows if row.get("run_id")})
    cohort_ids = sorted({row.get("cohort_id") for row in rows if row.get("cohort_id")})

    checksums = []
    for record in model_publications:
        source = args.model_reports / record["source_path"]
        if not source.is_file():
            raise SystemExit(f"model publication source was not staged: {source}")
        checksums.append((record["published_path"], str(source), digest(source)))
    args.checksums.write_text(
        "path\tsource_path\tsha256\n"
        + "".join(f"{published}\t{source}\t{checksum}\n" for published, source, checksum in sorted(checksums))
    )
    manifest = {
        "pipeline": "long-read",
        "configuration": config,
        "run_ids": run_ids,
        "cohort_ids": cohort_ids,
        "reference": {
            "path": str(args.reference_fasta),
            "sha256": digest(args.reference_fasta),
            "identity_sha256": digest(args.reference_identity),
        },
        "normalized_manifest": digest(args.normalized_manifest),
        "input_checksums": digest(args.input_checksums),
        "preparation_status": digest(args.preparation_status),
        "shard_plan_files": [{"path": str(path), "sha256": digest(path)} for path in sorted(Path("shard_plans").glob("*"))],
        "merge_provenance": merge_provenance,
        "stage_status": status,
        "final_product_checksum_file": digest(args.checksums),
    }
    args.run_manifest.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
