#!/usr/bin/env python3
"""Resolve ENA read-run evidence with a persistent, refreshable cache.

The cache stores the lossless response for each accession.  A failed optional
SRA audit is represented explicitly as ``unavailable`` in the normalized JSON.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path

FIELDS = ("run_accession,instrument_platform,instrument_model,fastq_ftp,fastq_md5,"
          "submitted_ftp,submitted_md5,submitted_format,library_strategy,library_source,library_selection")


def normalize_uri_list(value: str) -> str:
    """ENA often returns FTP paths without the ftp:// scheme."""
    return ";".join(
        ("ftp://" + item if item.startswith("ftp.") else item)
        for item in value.split(";") if item
    )


def ena_url(accession: str) -> str:
    query = urllib.parse.urlencode({"accession": accession, "result": "read_run", "fields": FIELDS, "format": "json"})
    return "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query


def ncbi_original_url(accession: str) -> str:
    return "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=" + urllib.parse.quote(accession)


def original_representation(filename: str) -> str:
    name = filename.lower()
    if name.endswith((".ccs.bam", ".hifi_reads.bam")): return "PACBIO_CCS_BAM"
    if name.endswith((".ccs.fastq", ".ccs.fq", ".ccs.fastq.gz", ".ccs.fq.gz")): return "PACBIO_CCS_FASTQ"
    if name.endswith(".subreads.bam"): return "PACBIO_SUBREAD_BAM"
    if any(marker in name for marker in (".flnc.", ".processed.", ".processed_reads.")): return "PACBIO_PROCESSED_FASTQ"
    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")): return "PACBIO_RAW_FASTQ"
    return "UNKNOWN"


def resolve_ncbi_original(accession: str, cache_dir: Path) -> dict:
    """Read NCBI's public Original Format inventory without downloading reads."""
    path = cache_dir / f"{accession}.ncbi.run.xml"
    try:
        if not path.exists():
            with urllib.request.urlopen(ncbi_original_url(accession), timeout=60) as response:
                partial = path.with_name(path.name + ".partial")
                partial.write_bytes(response.read())
                partial.replace(path)
        root = ET.fromstring(path.read_bytes())
        run = root.find(".//RUN")
        files = []
        for item in root.findall(".//SRAFile"):
            if item.attrib.get("supertype", "") != "Original":
                continue
            alternatives = [alt.attrib.get("url", "") for alt in item.findall("./Alternatives") if alt.attrib.get("url")]
            filename = item.attrib.get("filename", "")
            files.append({"filename": filename, "size": item.attrib.get("size", "0"),
                          "md5": item.attrib.get("md5", ""), "format": original_representation(filename),
                          "uri": alternatives[0] if alternatives else "", "supertype": "Original"})
        representations = [item["format"] for item in files if item["format"] != "UNKNOWN"]
        return {"ncbi_original_status": "available", "ncbi_original_alias": run.attrib.get("alias", "") if run is not None else "",
                "ncbi_original_files": files, "ncbi_original_representation": representations[0] if len(set(representations)) == 1 and representations else ("MIXED" if representations else "UNKNOWN"),
                "ncbi_original_xml": str(path)}
    except (OSError, ET.ParseError, urllib.error.URLError) as exc:
        return {"ncbi_original_status": "unavailable", "ncbi_original_alias": "", "ncbi_original_files": [],
                "ncbi_original_representation": "UNKNOWN", "ncbi_original_error": str(exc)[:240]}


def resolve(accession: str, cache_dir: Path, refresh: bool = False) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    raw_path = cache_dir / f"{accession}.ena.json"
    if raw_path.exists() and not refresh:
        raw = json.loads(raw_path.read_text())
    else:
        with urllib.request.urlopen(ena_url(accession), timeout=60) as response:
            raw_bytes = response.read()
        raw_path.write_bytes(raw_bytes)
        raw = json.loads(raw_bytes)
    row = raw[0] if isinstance(raw, list) and raw else raw
    if not row: raise ValueError(f"ENA returned no record for {accession}")
    audit = sra_audit(accession)
    for field in ("fastq_ftp", "submitted_ftp"):
        row[field] = normalize_uri_list(row.get(field, ""))
    original = resolve_ncbi_original(accession, cache_dir)
    return {"run_accession": accession, **row, **audit, **original, "ena_raw_response": str(raw_path)}


def sra_audit(accession: str) -> dict:
    command = shutil.which("sra-stat")
    if not command:
        return {"sra_platform": "unavailable", "sra_spot_group": "unavailable", "sra_status": "unavailable"}
    try:
        result = subprocess.run([command, accession], check=True, capture_output=True, text=True, timeout=120)
        text = result.stdout + result.stderr
        platform = next((line.split(":", 1)[1].strip() for line in text.splitlines() if "platform" in line.lower() and ":" in line), "unknown")
        spot = next((line.split(":", 1)[1].strip() for line in text.splitlines() if "spot group" in line.lower() and ":" in line), "unknown")
        return {"sra_platform": platform or "unknown", "sra_spot_group": spot or "unknown", "sra_status": "available"}
    except (OSError, subprocess.SubprocessError):
        return {"sra_platform": "unavailable", "sra_spot_group": "unavailable", "sra_status": "unavailable"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", help="one accession per line")
    parser.add_argument("cache_dir")
    parser.add_argument("output_json")
    parser.add_argument("--refresh_metadata", action="store_true")
    args = parser.parse_args()
    accessions = [line.strip() for line in Path(args.accessions).read_text().splitlines() if line.strip() and not line.startswith("#")]
    result = {accession: resolve(accession, Path(args.cache_dir), args.refresh_metadata) for accession in accessions}
    Path(args.output_json).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__": main()
