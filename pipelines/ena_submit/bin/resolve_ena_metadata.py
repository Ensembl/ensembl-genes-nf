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
ENA_SAMPLE_RE = re.compile(r"^(?:ERS|SRS|DRS)\d+$", re.IGNORECASE)
SAMPLE_RE = re.compile(r"^(?:SAMEA|SAMD|SAMN|ERS|SRS|DRS)\d+$", re.IGNORECASE)
EXPERIMENT_RE = re.compile(r"^(?:ERX|SRX|DRX)\d+$", re.IGNORECASE)


def split_values(value):
    return [x.strip() for x in re.split(r"[,;\s]+", value or "") if x.strip()]


def unique(values):
    return list(dict.fromkeys(x for x in values if x))


def accession_values(row, singular, plural, pattern):
    values = split_values(row.get(singular))
    values.extend(split_values(row.get(plural)))
    return unique(value for value in values if pattern.match(value))


def canonical_sample(record):
    """Return the ENA Sample accession used by SAMPLE_REF.

    ENA exposes two linked identifiers for a registered sample: an ENA Sample
    accession (ERS/SRS/DRS) and a BioSample accession (SAMEA/SAMN/SAMD). The
    analysis XML should reference the former; the latter is retained as
    provenance and can be used to inspect the parent BioSample record.
    """
    primary = (record.get("sample_accession") or "").strip()
    secondary = (record.get("secondary_sample_accession") or "").strip()
    if ENA_SAMPLE_RE.match(primary):
        return primary
    if ENA_SAMPLE_RE.match(secondary):
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
    metadata = [
        {field: (record.get(field) or "").strip() for field in FIELDS}
        for record in metadata
    ]

    source_sample_ids = unique(
        value
        for record in metadata
        for value in (record["sample_accession"], record["secondary_sample_accession"])
        if SAMPLE_RE.match(value)
    )
    expected_sample_ids = accession_values(
        row, "sample_accession", "sample_accessions", SAMPLE_RE
    )
    if expected_sample_ids and not set(expected_sample_ids).issubset(source_sample_ids):
        raise RuntimeError(
            "manifest sample accession(s) disagree with ENA: "
            f"manifest={expected_sample_ids}, ENA={source_sample_ids}"
        )

    source_experiment_ids = unique(
        record["experiment_accession"]
        for record in metadata
        if EXPERIMENT_RE.match(record["experiment_accession"])
    )
    expected_experiment_ids = accession_values(
        row, "experiment_accession", "experiment_accessions", EXPERIMENT_RE
    )
    if expected_experiment_ids and not set(expected_experiment_ids).issubset(source_experiment_ids):
        raise RuntimeError(
            "manifest experiment accession(s) disagree with ENA: "
            f"manifest={expected_experiment_ids}, ENA={source_experiment_ids}"
        )

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
    if row.get("sample_accession"):
        result["manifest_sample_accession"] = row["sample_accession"]
    if row.get("sample_accessions"):
        result["manifest_sample_accessions"] = row["sample_accessions"]
    if row.get("experiment_accession"):
        result["manifest_experiment_accession"] = row["experiment_accession"]
    if row.get("experiment_accessions"):
        result["manifest_experiment_accessions"] = row["experiment_accessions"]
    result["run_accessions"] = ",".join(runs)
    # These are authoritative ENA values.  Do not leave a valid-looking local
    # accession in the singular field when it disagrees with the API result.
    result["sample_accession"] = samples[0] if len(samples) == 1 else ""
    result["sample_accessions"] = ",".join(samples)
    result["biosample_accessions"] = ",".join(biosamples)
    result["secondary_sample_accessions"] = ",".join(
        unique(m["secondary_sample_accession"] for m in metadata)
    )
    result["experiment_accession"] = experiments[0] if len(experiments) == 1 else ""
    result["experiment_accessions"] = ",".join(experiments)
    result["source_study_accessions"] = ",".join(
        unique(m["study_accession"] for m in metadata)
    )
    result["secondary_study_accessions"] = ",".join(
        unique(m["secondary_study_accession"] for m in metadata)
    )
    result["scientific_name"] = species[0] if species else result.get("scientific_name", "")
    result["species"] = result.get("species") or (species[0].replace(" ", "_") if species else "")
    result["taxon_id"] = result.get("taxon_id") or (tax_ids[0] if tax_ids else "")
    for output_field, source_field in (
        ("library_strategies", "library_strategy"),
        ("library_sources", "library_source"),
        ("instrument_platforms", "instrument_platform"),
        ("instrument_models", "instrument_model"),
    ):
        result[output_field] = ",".join(unique(m[source_field] for m in metadata))
    # Use ENA's instrument platform for downstream analysis attributes when it
    # is available, while retaining the manifest value for audit.
    if row.get("platform"):
        result["manifest_platform"] = row["platform"]
    result["platform"] = result["instrument_platforms"] or result.get("platform", "")
    result["source_metadata_count"] = str(len(metadata))
    result["source_metadata_json"] = json.dumps(metadata, sort_keys=True, separators=(",", ":"))
    if len(samples) <= 1:
        result["aggregation_scope"] = "single_run" if len(runs) == 1 else "multi_run_single_sample"
    else:
        result["aggregation_scope"] = "multi_sample"
    result["metadata_source"] = "ena_portal_api"
    return result


def summarise_rows(rows):
    fields = {
        "run_accessions": [], "sample_accessions": [], "biosample_accessions": [],
        "secondary_sample_accessions": [],
        "experiment_accessions": [], "source_study_accessions": [],
        "secondary_study_accessions": [],
        "scientific_name": [], "taxon_id": [], "aggregation_scope": [],
        "library_strategies": [], "library_sources": [],
        "instrument_platforms": [], "instrument_models": [],
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
