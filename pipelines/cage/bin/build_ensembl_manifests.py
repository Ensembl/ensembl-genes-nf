#!/usr/bin/env python3
"""Build ENA CAGE sample/reference manifests using Ensembl species coverage."""

import argparse
import csv
import json
import re
import urllib.parse
import urllib.request
from urllib.error import HTTPError, URLError
from collections import Counter, defaultdict
from pathlib import Path


ENA_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "scientific_name",
    "tax_id", "library_name", "library_layout", "library_strategy",
    "library_source", "library_selection", "study_title", "first_public",
    "fastq_ftp", "submitted_ftp",
]


def fetch_json(url):
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ensembl-cage-manifest-builder/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def fetch_ena():
    params = {
        "result": "read_run",
        "query": 'library_selection="CAGE" AND library_source="TRANSCRIPTOMIC"',
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": "1000000",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "ensembl-cage-manifest-builder/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response:
        return list(csv.DictReader((line.decode("utf-8") for line in response), delimiter="\t"))


def normalise(value):
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def protocol_for(row):
    text = " ".join(row.get(field, "") for field in ("library_name", "study_title")).lower()
    if "heliscope" in text:
        return "helicoscope_cage"
    if "nano-cage" in text or "nanocage" in text:
        return "nano_cage"
    return "standard_cage"


def species_index():
    payloads = [fetch_json("https://rest.ensembl.org/info/species?content-type=application/json")]
    for division in ("EnsemblPlants", "EnsemblFungi", "EnsemblProtists", "EnsemblMetazoa", "EnsemblBacteria"):
        try:
            payloads.append(fetch_json(f"https://rest.ensembl.org/info/species?division={division};content-type=application/json"))
        except (HTTPError, URLError) as error:
            print(f"WARNING: could not query {division}: {error}")
    by_taxon = {}
    by_name = {}
    for payload in payloads:
      for species in payload["species"]:
        record = {
            "species": species.get("name", ""),
            "display_name": species.get("display_name", ""),
            "tax_id": str(species.get("taxon_id", "")),
            "assembly": species.get("assembly", ""),
            "accession": species.get("accession", ""),
            "division": species.get("division", ""),
        }
        if record["tax_id"]:
            by_taxon[record["tax_id"]] = record
        for name in [record["species"], record["display_name"]]:
            if name:
                by_name[normalise(name)] = record
    return by_taxon, by_name


def match_species(row, by_taxon, by_name):
    tax_id = str(row.get("tax_id", "")).strip()
    if tax_id in by_taxon:
        return by_taxon[tax_id], "tax_id"
    name = normalise(row.get("scientific_name", ""))
    if name in by_name:
        return by_name[name], "scientific_name"
    return None, "no_ensembl_match"


def normalise_layout(value):
    value = (value or "").upper()
    if value in {"SINGLE", "SE"}:
        return "SE"
    if value in {"PAIRED", "PE"}:
        return "PE"
    return "SE"


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--reference-root", required=True, help="Directory where Ensembl references will be staged later.")
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    by_taxon, by_name = species_index()
    ena_rows = fetch_ena()

    covered = []
    uncovered = []
    reference_rows = {}
    counts = Counter()
    for row in ena_rows:
        counts[row.get("scientific_name", "")] += 1
        reference, match_method = match_species(row, by_taxon, by_name)
        if reference is None or not reference["assembly"]:
            uncovered.append({**row, "coverage_status": "no_ensembl_assembly", "match_method": match_method})
            continue
        protocol = protocol_for(row)
        sample_id = row["run_accession"]
        fastq = row.get("fastq_ftp", "").split(";")[0]
        covered.append({
            "sample_id": sample_id,
            "species": reference["species"],
            "assembly": reference["assembly"],
            "assay_class": "cage",
            "protocol": protocol,
            "layout": normalise_layout(row.get("library_layout", "SE")),
            "accession": row["run_accession"],
            "fastq_1": "-",
            "fastq_2": "-",
            "tissue": "",
            "condition": "",
            "biological_replicate": "",
            "checksum_1": "-",
            "checksum_2": "-",
            "ena_species": row.get("scientific_name", ""),
            "study_accession": row.get("study_accession", ""),
            "first_public": row.get("first_public", ""),
            "fastq_ftp": row.get("fastq_ftp", ""),
            "match_method": match_method,
        })
        reference_rows[reference["species"]] = reference

    sample_fields = list(covered[0].keys()) if covered else ["sample_id"]
    write_tsv(output_dir / "cage_samples_ensembl.tsv", sample_fields, covered)
    write_tsv(output_dir / "cage_samples_uncovered.tsv", list(uncovered[0].keys()) if uncovered else ["run_accession"], uncovered)

    refs = []
    for species in sorted(reference_rows):
        reference = reference_rows[species]
        root = Path(args.reference_root) / reference["species"]
        refs.append({
            "assembly": reference["assembly"],
            "species": reference["species"],
            "fasta": str(root / "genome.fa"),
            "gtf": str(root / "annotation.gtf"),
            "chrom_sizes": str(root / "chrom.sizes"),
            "taxon_id": reference["tax_id"],
            "source": "Ensembl",
            "assembly_accession": reference["accession"],
            "reference_status": "needs_staging",
        })
    write_tsv(output_dir / "cage_references_ensembl.tsv", list(refs[0].keys()) if refs else ["assembly"], refs)

    references_by_species = {row["species"]: row for row in refs}
    samples_by_species = defaultdict(list)
    for row in covered:
        samples_by_species[row["species"]].append(row)
    species_index_rows = []
    for species in sorted(samples_by_species):
        species_dir = output_dir / "species" / species
        write_tsv(species_dir / "samples.tsv", sample_fields, samples_by_species[species])
        reference = references_by_species[species]
        write_tsv(species_dir / "reference.tsv", list(reference.keys()), [reference])
        species_index_rows.append({
            "species": species,
            "assembly": reference["assembly"],
            "sample_count": len(samples_by_species[species]),
            "samples_manifest": str(Path("species") / species / "samples.tsv"),
            "reference_manifest": str(Path("species") / species / "reference.tsv"),
        })
    write_tsv(output_dir / "species_manifest_index.tsv", list(species_index_rows[0].keys()) if species_index_rows else ["species"], species_index_rows)

    coverage = []
    all_species = sorted(set(counts) | {row["scientific_name"] for row in uncovered})
    covered_counts = Counter(row["ena_species"] for row in covered)
    for ena_species in all_species:
        matched = by_name.get(normalise(ena_species))
        coverage.append({
            "ena_species": ena_species,
            "run_count": counts[ena_species],
            "covered_run_count": covered_counts[ena_species],
            "status": "covered" if covered_counts[ena_species] else "not_covered",
            "ensembl_species": matched["species"] if matched else "",
            "ensembl_assembly": matched["assembly"] if matched else "",
            "reason": "" if covered_counts[ena_species] else "no Ensembl assembly match or assembly unavailable",
        })
    write_tsv(output_dir / "ensembl_coverage_report.tsv", list(coverage[0].keys()) if coverage else ["ena_species"], coverage)

    summary = {
        "ena_run_count": len(ena_rows),
        "covered_run_count": len(covered),
        "uncovered_run_count": len(uncovered),
        "covered_species_count": len(reference_rows),
        "survey_species_count": len(counts),
    }
    (output_dir / "manifest_build_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
