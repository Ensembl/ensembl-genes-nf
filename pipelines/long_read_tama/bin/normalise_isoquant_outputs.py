#!/usr/bin/env python3
"""Copy version-dependent IsoQuant products to stable pipeline names."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def one(root: Path, suffix: str) -> Path:
    matches = sorted(path for path in root.rglob(f"*{suffix}") if path.is_file())
    if not matches:
        raise SystemExit(f"IsoQuant output is missing required product '*{suffix}': {root}")
    if len(matches) > 1:
        raise SystemExit(f"IsoQuant output has multiple products matching '*{suffix}': {matches}")
    return matches[0]


def optional(root: Path, suffix: str) -> Path | None:
    matches = sorted(path for path in root.rglob(f"*{suffix}") if path.is_file())
    if len(matches) > 1:
        raise SystemExit(f"IsoQuant output has multiple products matching '*{suffix}': {matches}")
    return matches[0] if matches else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("destination")
    parser.add_argument("scope")
    args = parser.parse_args()

    source = Path(args.source)
    destination = Path(args.destination)
    destination.mkdir(parents=True, exist_ok=False)

    model_reads = one(source, ".transcript_model_reads.tsv.gz")
    read_info = optional(source, ".read_info.tsv.gz")
    products = {
        "transcript_models.gtf": one(source, ".transcript_models.gtf"),
        "transcript_model_reads.tsv.gz": model_reads,
    }
    read_info_source = "native_read_info"
    if read_info:
        products["read_info.tsv.gz"] = read_info
    else:
        # IsoQuant's annotation-free mode does not create read_info.tsv.gz;
        # transcript_model_reads is the native per-read evidence in that mode.
        products["read_info.tsv.gz"] = model_reads
        read_info_source = "transcript_model_reads_fallback"
    for suffix in (
        ".discovered_transcript_counts.tsv",
        ".discovered_transcript_tpm.tsv",
        ".discovered_gene_counts.tsv",
        ".discovered_gene_tpm.tsv",
    ):
        match = optional(source, suffix)
        if match:
            products[suffix.lstrip(".")] = match

    rows = []
    for stable_name, source_path in products.items():
        target = destination / stable_name
        shutil.copy2(source_path, target)
        rows.append((stable_name, str(source_path.relative_to(source))))

    with (destination / "product_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["scope", "stable_name", "source_name", "source_type"])
        writer.writerows(
            (args.scope, stable, source_name, read_info_source if stable == "read_info.tsv.gz" else "native")
            for stable, source_name in rows
        )


if __name__ == "__main__":
    main()
