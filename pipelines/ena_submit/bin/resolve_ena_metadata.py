#!/usr/bin/env python3
"""Resolve and validate source metadata for ENA alignment submissions."""

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession",
    "experiment_accession", "study_accession", "secondary_study_accession",
    "scientific_name", "tax_id", "library_strategy", "library_source",
    "instrument_platform", "instrument_model",
]
RUN_RE = re.compile(r"^(ERR|SRR|DRR)\d+$", re.IGNORECASE)
BIO_SAMPLE_RE = re.compile(r"^SAM(?:EA|D|N)\d+$", re.IGNORECASE)


def split_values(value):
    return [x.strip() for x in re.split(r"[,;\s]+", value or "") if x.strip()]


def unique(values):
    return list(dict.fromkeys(x for x in values if x))


def canonical_sample(record):
    """Prefer the EMBL-EBI BioSample accession over an ENA sample alias."""
    primary = (record.get("sample_accession") or "").strip()
    secondary = (record.get("secondary_sample_accession") or "").strip()
    if BIO_SAMPLE_RE.match(primary):
        return primary
    if BIO_SAMPLE_RE.match(secondary):
        return secondary
    return primary or secondary


def fetch_run(accession, api_url):
    query = urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
    })
    request = Request(f"{api_url}?{query}", headers={"Accept": "text/tab-separated-values"})
    with urlopen(request, timeout=30) as response:
        rows = list(csv.DictReader((line.decode("utf-8") for line in response), delimiter="\t"))
    if len(rows) != 1:
        raise RuntimeError(f"ENA returned {len(rows)} records for {accession}")
    return rows[0]


def read_run_list(row):
    runs = split_values(row.get("run_accessions") or row.get("run_accession"))
    path = (row.get("run_list_path") or "").strip()
    if path:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for source in reader:
                runs.extend(split_values(source.get("run_accession") or source.get("run")))
    return unique(runs)


def resolve_row(row, api_url, strict=True):
    runs = read_run_list(row)
    invalid = [run for run in runs if not RUN_RE.match(run)]
    if invalid:
        raise RuntimeError(f"invalid ENA run accession(s): {', '.join(invalid)}")
    if not runs:
        if strict:
            raise RuntimeError("row has no source run accession")
        return row

    metadata = [fetch_run(run, api_url) for run in runs]
    samples = unique([canonical_sample(m) for m in metadata])
    biosamples = unique([
        value
        for record in metadata
        for value in (record.get("sample_accession", ""), record.get("secondary_sample_accession", ""))
        if BIO_SAMPLE_RE.match((value or "").strip())
    ])
    experiments = unique([m.get("experiment_accession", "") for m in metadata])
    species = unique([m.get("scientific_name", "") for m in metadata])
    tax_ids = unique([m.get("tax_id", "") for m in metadata])
    if len(species) > 1 or len(tax_ids) > 1:
        raise RuntimeError(f"source runs disagree on taxonomy: species={species}, tax_id={tax_ids}")

    requested_taxon = (row.get("taxon_id") or "").strip()
    if requested_taxon and tax_ids and requested_taxon != tax_ids[0]:
        raise RuntimeError(f"manifest taxon_id {requested_taxon} disagrees with ENA tax_id {tax_ids[0]}")
    requested_species = (row.get("species") or "").strip()
    if requested_species and species and requested_species.lower() not in {species[0].lower(), species[0].replace(" ", "_").lower()}:
        raise RuntimeError(f"manifest species {requested_species} disagrees with ENA scientific_name {species[0]}")

    result = dict(row)
    result["run_accessions"] = ",".join(runs)
    result["sample_accessions"] = ",".join(samples)
    result["biosample_accessions"] = ",".join(biosamples)
    result["experiment_accessions"] = ",".join(experiments)
    result["source_study_accessions"] = ",".join(unique([m.get("study_accession", "") for m in metadata]))
    result["scientific_name"] = species[0] if species else result.get("scientific_name", "")
    result["species"] = result.get("species") or (species[0].replace(" ", "_") if species else "")
    result["taxon_id"] = result.get("taxon_id") or (tax_ids[0] if tax_ids else "")
    if len(samples) <= 1:
        result["aggregation_scope"] = "single_run" if len(runs) == 1 else "multi_run_single_sample"
    else:
        result["aggregation_scope"] = "multi_sample"
    result["metadata_source"] = "ena_portal_api"
    return result


def summarise_rows(rows):
    fields = {
        "run_accessions": [], "sample_accessions": [], "biosample_accessions": [],
        "experiment_accessions": [], "source_study_accessions": [],
        "scientific_name": [], "taxon_id": [], "aggregation_scope": [],
    }
    for row in rows:
        for field in fields:
            if field == "scientific_name":
                value = (row.get(field) or "").strip()
                if value:
                    fields[field].append(value)
            else:
                fields[field].extend(split_values(row.get(field)))
    summary = {field: ",".join(unique(values)) for field, values in fields.items()}
    summary["metadata_source"] = "ena_portal_api"
    summary["species"] = summary["scientific_name"].replace(" ", "_")
    if rows:
        for field in ("project_alias", "assembly_accession", "release", "partial_release_label", "last_geneset_update"):
            if rows[0].get(field):
                summary[field] = rows[0][field]
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary-output")
    parser.add_argument("--project-alias", default="")
    parser.add_argument("--analysis-id", default="")
    parser.add_argument("--assembly-accession", default="")
    parser.add_argument("--release", default="")
    parser.add_argument("--api-url", default="https://www.ebi.ac.uk/ena/portal/api/filereport")
    parser.add_argument("--allow-missing-runs", action="store_true")
    args = parser.parse_args()

    with open(args.input, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError("metadata input is empty")

    resolved = [resolve_row(row, args.api_url, strict=not args.allow_missing_runs) for row in rows]
    columns = list(dict.fromkeys(key for row in resolved for key in row))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(resolved)
    if args.summary_output:
        summary = summarise_rows(resolved)
        summary["project_alias"] = args.project_alias or summary.get("project_alias", "")
        summary["analysis_id"] = args.analysis_id
        summary["assembly_accession"] = args.assembly_accession or summary.get("assembly_accession", "")
        summary["release"] = args.release or summary.get("release", "")
        with open(args.summary_output, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t")
            writer.writeheader()
            writer.writerow(summary)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
