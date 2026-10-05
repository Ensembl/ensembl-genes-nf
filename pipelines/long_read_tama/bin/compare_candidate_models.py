#!/usr/bin/env python3
"""Create a backend-neutral structural comparison for canonical BED12 files.

The report deliberately compares models without merging them.  Native model
products and read-level evidence remain the source of truth for each backend;
this script only computes comparable structural summaries and pairwise overlap.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Model:
    name: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]

    @property
    def intron_chain(self) -> tuple[tuple[int, int], ...]:
        return tuple((self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1))


def parse_bed12(path: Path) -> list[Model]:
    models: list[Model] = []
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith("track") or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 12:
            raise SystemExit(f"{path}:{line_number}: expected BED12, found {len(fields)} columns")
        chrom, start_text, end_text, name, _score, strand = fields[:6]
        block_count = int(fields[9])
        block_sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
        block_starts = [int(value) for value in fields[11].rstrip(",").split(",") if value]
        if block_count != len(block_sizes) or block_count != len(block_starts):
            raise SystemExit(f"{path}:{line_number}: BED12 block count does not match block lists")
        start = int(start_text)
        end = int(end_text)
        exons = tuple((start + offset, start + offset + size) for offset, size in zip(block_starts, block_sizes))
        if not exons or exons[0][0] < start or exons[-1][1] > end or any(a >= b for a, b in exons):
            raise SystemExit(f"{path}:{line_number}: invalid BED12 exon coordinates")
        models.append(Model(name=name, chrom=chrom, strand=strand, exons=exons))
    if not models:
        raise SystemExit(f"{path}: no BED12 models")
    return models


def backend_name(path: Path) -> str:
    suffix = "_combined_models.bed"
    if path.name.endswith(suffix):
        return path.name[: -len(suffix)]
    return path.stem


def summarise(models: list[Model]) -> dict[str, object]:
    chains = {(model.chrom, model.strand, model.intron_chain) for model in models}
    single_exon = sum(len(model.exons) == 1 for model in models)
    lengths = [end - start for model in models for start, end in [(model.exons[0][0], model.exons[-1][1])]]
    return {
        "model_count": len(models),
        "unique_intron_chain_count": len(chains),
        "duplicate_model_count": len(models) - len(chains),
        "single_exon_count": single_exon,
        "multi_exon_count": len(models) - single_exon,
        "chromosome_count": len({model.chrom for model in models}),
        "mean_span": sum(lengths) / len(lengths),
        "median_span": sorted(lengths)[len(lengths) // 2],
        "min_span": min(lengths),
        "max_span": max(lengths),
    }


def structural_keys(models: list[Model]) -> set[tuple[str, str, tuple[tuple[int, int], ...]]]:
    return {(model.chrom, model.strand, model.intron_chain) for model in models}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("models", nargs="+", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--stage", default="unknown")
    args = parser.parse_args()

    parsed = {backend_name(path): parse_bed12(path) for path in args.models}
    if len(parsed) != len(args.models):
        raise SystemExit("model filenames do not identify unique backends")
    keys = {backend: structural_keys(models) for backend, models in parsed.items()}
    rows: list[dict[str, object]] = []
    for backend, models in sorted(parsed.items()):
        digest = hashlib.sha256(next(path for path in args.models if backend_name(path) == backend).read_bytes()).hexdigest()
        rows.append({"backend": backend, "model_sha256": digest, **summarise(models)})

    pairwise: list[dict[str, object]] = []
    backends = sorted(keys)
    for index, left in enumerate(backends):
        for right in backends[index + 1 :]:
            intersection = len(keys[left] & keys[right])
            union = len(keys[left] | keys[right])
            pairwise.append({
                "left_backend": left,
                "right_backend": right,
                "shared_unique_intron_chains": intersection,
                "union_unique_intron_chains": union,
                "jaccard_unique_intron_chains": intersection / union if union else 0.0,
            })

    args.output_tsv.write_text(
        "backend\tmodel_sha256\tmodel_count\tunique_intron_chain_count\tduplicate_model_count\t"
        "single_exon_count\tmulti_exon_count\tchromosome_count\tmean_span\tmedian_span\tmin_span\tmax_span\n"
        + "\n".join(
            "\t".join(str(row[key]) for key in (
                "backend", "model_sha256", "model_count", "unique_intron_chain_count", "duplicate_model_count",
                "single_exon_count", "multi_exon_count", "chromosome_count", "mean_span", "median_span", "min_span", "max_span"
            )) for row in rows
        ) + "\n"
    )
    args.output_json.write_text(json.dumps({"backends": rows, "pairwise": pairwise}, indent=2) + "\n")
    if args.manifest:
        manifest_rows = []
        for backend, models in sorted(parsed.items()):
            bed_path = next(path for path in args.models if backend_name(path) == backend)
            digest = hashlib.sha256(bed_path.read_bytes()).hexdigest()
            for model in models:
                introns = ";".join(f"{left}-{right}" for left, right in model.intron_chain) or "NONE"
                manifest_rows.append({
                    "backend": backend,
                    "stage": args.stage,
                    "canonical_model_id": model.name,
                    "native_model_id": model.name,
                    "chromosome": model.chrom,
                    "strand": model.strand,
                    "start": model.exons[0][0],
                    "end": model.exons[-1][1],
                    "exon_count": len(model.exons),
                    "intron_chain": introns,
                    "canonical_bed": bed_path.name,
                    "canonical_bed_sha256": digest,
                    "read_support": "UNKNOWN",
                    "accession_support": "UNKNOWN",
                    "provenance_status": "STRUCTURAL_ONLY",
                })
        fields = [
            "backend", "stage", "canonical_model_id", "native_model_id", "chromosome", "strand",
            "start", "end", "exon_count", "intron_chain", "canonical_bed", "canonical_bed_sha256",
            "read_support", "accession_support", "provenance_status",
        ]
        with args.manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
