#!/usr/bin/env python3
"""Discover and select transcriptomic long-read evidence without bulk download."""
from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from read_input_classification import probe_fastq
from resolve_metadata import resolve_ncbi_original

FIELDS = ",".join([
    "run_accession", "sample_accession", "experiment_accession", "study_accession",
    "instrument_platform", "instrument_model", "library_strategy", "library_source",
    "library_selection", "library_layout", "read_count", "base_count", "fastq_bytes", "fastq_ftp",
    "fastq_md5", "submitted_ftp", "submitted_md5", "submitted_format",
])
BIO_FIELDS = ("tissue", "organism_part", "cell_type", "cell_line", "developmental_stage",
              "life_stage", "age", "sex", "strain", "breed", "isolate", "description",
              "sample_title", "name")


def get_json(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": "ensembl-long-read-tama/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def cached_json(url: str, path: Path):
    if path.exists():
        return json.loads(path.read_text())
    data = get_json(url)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
    return data


def first_characteristic(data: dict, keys: tuple[str, ...]) -> str:
    chars = data.get("characteristics", {}) if isinstance(data, dict) else {}
    for key in keys:
        values = chars.get(key) or []
        if values:
            value = values[0].get("text", values[0]) if isinstance(values[0], dict) else values[0]
            if value:
                return str(value)
    return "unknown"


def normalise_tissue(raw: str) -> str:
    text = raw.lower().replace("_", " ").strip()
    if text == "unknown": return text
    if "," in text or "mixed" in text or "whole body" in text: return "mixed_tissues"
    aliases = {"heart": "heart", "liver": "liver", "brain": "brain", "testis": "testis", "testes": "testis",
               "ovary": "ovary", "ovaries": "ovary", "muscle": "muscle", "kidney": "kidney",
               "lung": "lung", "blood": "blood", "skin": "skin", "embryo": "embryo"}
    for key, value in aliases.items():
        if key in text: return value
    return "other:" + re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def artifact(urls: str, md5s: str, fmt: str = "") -> tuple[str, str, str]:
    entries = [x for x in urls.split(";") if x]
    checksums = [x for x in md5s.split(";") if x]
    for uri, checksum in zip(entries, checksums):
        name = Path(urllib.parse.urlparse(uri if "://" in uri else "ftp://" + uri).path).name
        lower = name.lower()
        if lower.endswith(".subreads.bam"): return "PACBIO_SUBREAD_BAM", uri, checksum
        if lower.endswith((".ccs.bam", ".hifi_reads.bam")): return "PACBIO_CCS_BAM", uri, checksum
    return "", "", ""


def _first_number(value: str) -> int:
    try:
        return int((value or "0").split(";")[0])
    except (TypeError, ValueError):
        return 0


def _representation_penalty(representation: str) -> int:
    return {
        "ONT_FASTQ": 0,
        "PACBIO_CCS_BAM": 0,
        "PACBIO_CCS_FASTQ": 1,
        "PACBIO_PROCESSED_BAM": 1,
        "PACBIO_PROCESSED_FASTQ": 2,
        "PACBIO_SUBREAD_BAM": 2,
        "PACBIO_RAW_SUBREAD_FASTQ_GROUPABLE": 3,
        "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE": 5,
        "PACBIO_UNKNOWN": 4,
    }.get(representation, 4)


def _technical_annotation(row: dict) -> tuple[int, str, str]:
    representation = row["representation_class"]
    penalty = _representation_penalty(representation)
    if "UNGROUPABLE" in representation:
        label = "high"
        warning = "raw PacBio FASTQ lacks observable molecule grouping; subreads must not be treated as independent molecules"
    elif representation in {"PACBIO_UNKNOWN", "PROBE_FAILED"}:
        label = "high"
        warning = "PacBio representation remains unresolved; reviewer must confirm artifact and molecule identity"
    elif "RAW_SUBREAD" in representation:
        label = "medium-high"
        warning = "raw PacBio subreads require molecule-aware processing before annotation evidence is counted"
    else:
        label = "low" if penalty == 0 else "medium"
        warning = ""
    return penalty, label, warning


def discover(taxon: str, tree: bool, cache: Path, probe: bool, target: int,
             soft_download_budget: float = 250.0, soft_raw_subread_budget: int = 1):
    query_taxon = f"tax_tree({taxon})" if tree else f"tax_eq({taxon})"
    query = f"{query_taxon} AND library_source=TRANSCRIPTOMIC AND (instrument_platform=OXFORD_NANOPORE OR instrument_platform=PACBIO_SMRT)"
    params = urllib.parse.urlencode({"display": "report", "domain": "read", "result": "read_run", "query": query, "fields": FIELDS})
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + params
    raw = cached_json(url, cache / "ena_discovery.json") if (cache / "ena_discovery.json").exists() else None
    if raw is None:
        request = urllib.request.Request(url, headers={"User-Agent": "ensembl-long-read-tama/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            lines = response.read().decode().splitlines()
        reader = csv.DictReader(lines, delimiter="\t")
        raw = list(reader)
        (cache / "ena_discovery.json").parent.mkdir(parents=True, exist_ok=True)
        (cache / "ena_discovery.json").write_text(json.dumps(raw, indent=2) + "\n")
    rows = []
    for row in raw:
        sample = (row.get("sample_accession") or "unknown").split(";")[0]
        bio_url = "https://www.ebi.ac.uk/biosamples/samples/" + urllib.parse.quote(sample)
        bio = cached_json(bio_url, cache / "biosamples" / f"{sample}.json") if sample != "unknown" else {}
        raw_tissue = first_characteristic(bio, ("tissue", "organism_part"))
        subread_class, subread_uri, subread_md5 = artifact(row.get("submitted_ftp", ""), row.get("submitted_md5", ""), row.get("submitted_format", ""))
        original = resolve_ncbi_original(row.get("run_accession", ""), cache / "ncbi_original")
        original_representation = original.get("ncbi_original_representation", "UNKNOWN")
        observed = "NOT_PROBED"
        if row.get("instrument_platform") == "OXFORD_NANOPORE": observed = "ONT"
        elif subread_class: observed = subread_class
        elif original_representation in {"PACBIO_CCS_FASTQ", "PACBIO_CCS_BAM", "PACBIO_SUBREAD_BAM", "PACBIO_PROCESSED_FASTQ"}:
            observed = original_representation
        elif probe and row.get("fastq_ftp"):
            fastq = row["fastq_ftp"].split(";")[0]
            if fastq.startswith("ftp."): fastq = "ftp://" + fastq
            try: observed = probe_fastq(fastq, limit=100).get("header_representation", "UNKNOWN")
            except (OSError, ValueError, TimeoutError) as exc: observed = "PROBE_FAILED:" + str(exc)[:120]
        fastq_name = Path(urllib.parse.urlparse(row.get("fastq_ftp", "").split(";")[0]).path).name.lower()
        if original_representation in {"PACBIO_CCS_FASTQ", "PACBIO_CCS_BAM", "PACBIO_SUBREAD_BAM", "PACBIO_PROCESSED_FASTQ"}: representation = original_representation
        elif observed == "PACBIO_SUBREAD": representation = "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE"
        elif observed == "PACBIO_CCS": representation = "PACBIO_CCS_FASTQ"
        elif observed == "ONT": representation = "ONT_FASTQ"
        elif subread_class: representation = subread_class
        elif row.get("instrument_platform") == "PACBIO_SMRT" and "subread" in fastq_name: representation = "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE"
        else: representation = "PACBIO_UNKNOWN" if row.get("instrument_platform") == "PACBIO_SMRT" else "ONT_FASTQ"
        selected_uri = subread_uri or row.get("fastq_ftp", "").split(";")[0]
        selected_md5 = subread_md5 or row.get("fastq_md5", "").split(";")[0]
        selected_name = Path(urllib.parse.urlparse(selected_uri if "://" in selected_uri else "ftp://" + selected_uri).path).name
        rows.append({"species_taxon_id": taxon, "sample_accession": sample, "run_accession": row.get("run_accession", ""),
                     "experiment_accession": row.get("experiment_accession", "unknown"), "study_accession": row.get("study_accession", "unknown"),
                     "raw_tissue": raw_tissue, "normalised_tissue": normalise_tissue(raw_tissue),
                     "developmental_stage": first_characteristic(bio, ("developmental_stage", "life_stage", "age")),
                     "sex": first_characteristic(bio, ("sex",)), "strain_or_breed": first_characteristic(bio, ("strain", "breed", "isolate")),
                     "description": first_characteristic(bio, ("description", "sample_title", "name")),
                     "platform": row.get("instrument_platform", "unknown"), "instrument_model": row.get("instrument_model", "unknown"),
                     "library_strategy": row.get("library_strategy", "unknown"), "library_source": row.get("library_source", "unknown"),
                     "read_count": row.get("read_count", "0"), "base_count": row.get("base_count", "0"),
                     "fastq_bytes": row.get("fastq_bytes", "0"),
                     "fastq_ftp": row.get("fastq_ftp", ""), "fastq_md5": row.get("fastq_md5", ""),
                     "submitted_ftp": row.get("submitted_ftp", ""), "submitted_md5": row.get("submitted_md5", ""),
                     "submitted_format": row.get("submitted_format", ""), "header_representation": observed,
                     "ncbi_original_alias": original.get("ncbi_original_alias", ""), "ncbi_original_representation": original_representation,
                     "ncbi_original_files": json.dumps(original.get("ncbi_original_files", []), separators=(",", ":")),
                     "representation_class": representation, "selected_artifact": selected_uri,
                     "selected_artifact_md5": selected_md5, "tissue": raw_tissue, "filename": selected_name,
                     "url": selected_uri, "md5": selected_md5, "source": "ENA"})
    return select(rows, target, soft_download_budget, soft_raw_subread_budget)


def select(rows: list[dict], target: int, soft_download_budget: float, soft_raw_subread_budget: int):
    by_sample = defaultdict(list)
    for row in rows: by_sample[row["sample_accession"]].append(row)
    # Select one technically preferred representative per biological sample,
    # while retaining every run in the full inventory for review.
    candidates = []
    for group in by_sample.values():
        candidates.append(max(group, key=lambda r: (
            -_representation_penalty(r["representation_class"]),
            _first_number(r.get("base_count", "0")),
        )))
    selected, represented = [], set()
    represented_stages, represented_contexts = set(), set()
    used_download_gb, used_raw_subreads = 0.0, 0
    while candidates and len(selected) < target:
        def score(row):
            tissue = row["normalised_tissue"]
            stage = row["developmental_stage"]
            context = (tissue, stage)
            biological = (100 if tissue not in represented else 0)
            biological += (55 if stage != "unknown" and stage not in represented_stages else 0)
            biological += (25 if context not in represented_contexts else 0)
            if tissue == "mixed_tissues":
                biological += 20
            evidence = min(25, max(0, _first_number(row.get("read_count", "0")) // 100000))
            technical_level, _, _ = _technical_annotation(row)
            technical = technical_level * 16
            size_gb = _first_number(row.get("fastq_bytes", "0")) / 1e9
            acquisition = min(35, size_gb / 10.0)
            raw_over_budget = max(0, used_raw_subreads + ("RAW_SUBREAD" in row["representation_class"]) - soft_raw_subread_budget)
            budget_over = max(0.0, used_download_gb + size_gb - soft_download_budget)
            return biological + evidence - technical - acquisition - raw_over_budget * 35 - budget_over / 5
        row = max(candidates, key=score); candidates.remove(row)
        gain = (100 if row["normalised_tissue"] not in represented else 0)
        gain += 55 if row["developmental_stage"] != "unknown" and row["developmental_stage"] not in represented_stages else 0
        gain += 25 if (row["normalised_tissue"], row["developmental_stage"]) not in represented_contexts else 0
        penalty, penalty_label, warning = _technical_annotation(row)
        size_gb = _first_number(row.get("fastq_bytes", "0")) / 1e9
        row["technical_penalty"] = penalty_label
        row["biological_gain"] = gain
        row["warnings"] = warning
        row["selection_status"] = "SELECTED_WITH_WARNING" if warning or used_download_gb + size_gb > soft_download_budget or used_raw_subreads + ("RAW_SUBREAD" in row["representation_class"]) > soft_raw_subread_budget else "SELECTED"
        row["selection_reason"] = (f"Adds {row['normalised_tissue']} / {row['developmental_stage']} evidence; "
                                    f"representation {row['representation_class']}, estimated FASTQ {size_gb:.1f} GB, "
                                    f"technical penalty {penalty_label}")
        selected.append(row); represented.add(row["normalised_tissue"])
        represented_stages.add(row["developmental_stage"])
        represented_contexts.add((row["normalised_tissue"], row["developmental_stage"]))
        used_download_gb += size_gb
        used_raw_subreads += int("RAW_SUBREAD" in row["representation_class"])
    for row in rows:
        if "technical_penalty" not in row:
            _, row["technical_penalty"], row["warnings"] = _technical_annotation(row)
        if "selection_status" not in row:
            row["selection_status"] = "NOT_SELECTED"; row["biological_gain"] = 0
            row["selection_reason"] = "Biological context already represented or lower marginal value after sample-level comparison"
    return rows, selected


def main():
    p = argparse.ArgumentParser(); p.add_argument("taxon_id"); p.add_argument("output_dir"); p.add_argument("--tree", action="store_true"); p.add_argument("--cache-dir", default=None); p.add_argument("--no-probe", action="store_true"); p.add_argument("--target", type=int, default=8); p.add_argument("--soft-download-budget", type=float, default=250.0); p.add_argument("--soft-raw-subread-budget", type=int, default=1)
    a = p.parse_args(); out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True); cache = Path(a.cache_dir or out / "cache")
    rows, selected = discover(a.taxon_id, a.tree, cache, not a.no_probe, a.target, a.soft_download_budget, a.soft_raw_subread_budget)
    fields = list(rows[0]) if rows else ["run_accession", "selection_status"]
    with (out / "long_read_inventory.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(rows)
    with (out / "proposed_long_read_manifest.tsv").open("w", newline="") as h:
        # This is the versioned manifest shape accepted by normalise_manifest.py.
        # It contains selected runs only; selection metadata remain in inventory.
        fields2 = ["run_accession", "tissue", "filename", "url", "md5", "source", "platform", "description", "selection_status", "selection_reason"]
        w = csv.DictWriter(h, fieldnames=fields2, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(selected)
    summary = {"discovered_runs": len(rows), "biological_samples": len({r['sample_accession'] for r in rows}),
               "selected_runs": len(selected), "soft_download_budget_gb": a.soft_download_budget,
               "soft_raw_subread_budget": a.soft_raw_subread_budget,
               "selected_with_warning": sum(r["selection_status"] == "SELECTED_WITH_WARNING" for r in selected)}
    (out / "selection_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__": main()
