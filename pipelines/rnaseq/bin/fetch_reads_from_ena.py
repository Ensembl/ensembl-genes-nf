#!/usr/bin/env python3
"""
Fetch short-read RNA-seq data from ENA by BioProject or run accessions.

Queries the ENA Portal API to resolve BioProject → run accessions,
downloads FASTQ files via ENA FTP, and writes a sample_sheet.csv
compatible with the rnaseq pipeline.

Usage:
  fetch_reads_from_ena.py --accession PRJEB12345 --outdir .
  fetch_reads_from_ena.py --accession ERR123456,ERR123457 --outdir .
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import requests

ENA_PORTAL_URL  = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_SEARCH_URL  = "https://www.ebi.ac.uk/ena/portal/api/search"
BIOSAMPLE_URL   = "https://www.ebi.ac.uk/biosamples/samples"
ENA_FTP_BASE    = "ftp://ftp.sra.ebi.ac.uk/vol1"

MAX_RETRIES = 3
RETRY_DELAY = 10  # seconds


def _get(url: str, params: dict, retries: int = MAX_RETRIES) -> dict | list:
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=120)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:
            if attempt == retries - 1:
                raise RuntimeError(f"ENA API request failed after {retries} attempts: {exc}") from exc
            print(f"  Retry {attempt+1}/{retries} for {url}: {exc}", file=sys.stderr)
            time.sleep(RETRY_DELAY)


def resolve_runs(accession: str) -> list[dict]:
    """
    Resolve a BioProject or comma-separated run accessions to a list of run dicts.
    Each dict has: run_accession, fastq_ftp, library_layout, sample_accession.
    """
    # Comma-separated run accessions (SRR/ERR/DRR)
    if re.match(r'^[SED]RR\d+', accession):
        runs_str = accession
        query = f"run_accession={runs_str.replace(',', '%20OR%20run_accession=')}"
        params = {
            "accession": accession.split(',')[0],  # just use first for API call
            "result":    "read_run",
            "fields":    "run_accession,fastq_ftp,library_layout,sample_accession,instrument_platform",
            "format":    "json",
            "limit":     500,
        }
    elif re.match(r'^PRJ[NE][A-Z]\d+', accession) or re.match(r'^ERP\d+|SRP\d+', accession):
        params = {
            "accession": accession,
            "result":    "read_run",
            "fields":    "run_accession,fastq_ftp,library_layout,sample_accession,instrument_platform",
            "format":    "json",
            "limit":     1000,
        }
    else:
        raise ValueError(f"Unrecognised accession format: {accession!r}")

    data = _get(ENA_PORTAL_URL, params)
    if not isinstance(data, list):
        raise RuntimeError(f"Unexpected ENA API response: {data}")
    return data


def get_tissue_label(sample_accession: str) -> str:
    """Query BioSamples for tissue/organ label to use as the sample id."""
    try:
        resp = requests.get(f"{BIOSAMPLE_URL}/{sample_accession}", timeout=30)
        if not resp.ok:
            return sample_accession
        info = resp.json()
        chars = info.get("characteristics", {})
        for key in ("tissue", "organism part", "cell_type", "source name"):
            if key in chars:
                label = chars[key][0]["text"].lower()
                label = re.sub(r"[ ;\(\)\/\\]", "_", label)
                label = re.sub(r"[^a-z0-9_]", "", label)
                return label or sample_accession
        return sample_accession
    except Exception:
        return sample_accession


def download_fastq(ftp_path: str, outdir: Path, retries: int = MAX_RETRIES) -> Optional[Path]:
    """Download one FASTQ file from ENA FTP using wget. Returns local path or None on failure."""
    # ftp_path may be like: ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR123456/ERR123456_1.fastq.gz
    url = f"ftp://{ftp_path}" if not ftp_path.startswith("ftp://") else ftp_path
    filename = url.split("/")[-1]
    local = outdir / filename
    if local.exists():
        print(f"  Already downloaded: {filename}", file=sys.stderr)
        return local
    cmd = ["wget", "-q", f"--tries={retries}", "--timeout=600", url, "-O", str(local)]
    print(f"  Downloading: {filename}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        print(f"  WARN: download failed for {url}", file=sys.stderr)
        if local.exists():
            local.unlink()
        return None
    return local


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession",    required=True,  help="BioProject or comma-separated run accessions")
    parser.add_argument("--outdir",       required=True,  help="Output directory for FASTQ files")
    parser.add_argument("--max-runs",     type=int, default=50, help="Maximum number of runs to download")
    parser.add_argument("--strandedness", default="auto",  help="Strandedness override (forward/reverse/unstranded/auto)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve accession → run list
    print(f"Resolving ENA accession: {args.accession}", file=sys.stderr)
    runs = resolve_runs(args.accession)

    # Filter to Illumina short-read only
    runs = [r for r in runs if r.get("instrument_platform", "ILLUMINA") == "ILLUMINA"]

    if not runs:
        print("ERROR: no Illumina short-read runs found for this accession", file=sys.stderr)
        sys.exit(1)

    if len(runs) > args.max_runs:
        print(f"  Capping at {args.max_runs} of {len(runs)} runs (--max-runs)", file=sys.stderr)
        runs = runs[:args.max_runs]

    rows = []
    for run in runs:
        acc    = run["run_accession"]
        layout = run.get("library_layout", "PAIRED").upper()
        ftp    = run.get("fastq_ftp", "")
        if not ftp:
            print(f"  WARN: no FTP URL for {acc}, skipping", file=sys.stderr)
            continue

        ftp_files = [f.strip() for f in ftp.split(";") if f.strip()]

        # Get tissue label for id
        sample_id = run.get("sample_accession", acc)
        label     = get_tissue_label(sample_id)
        safe_id   = f"{label}_{acc}"

        # Download
        local_files = []
        for f in ftp_files:
            lf = download_fastq(f, outdir)
            if lf:
                local_files.append(str(lf))

        if not local_files:
            print(f"  WARN: no files downloaded for {acc}, skipping", file=sys.stderr)
            continue

        strand = args.strandedness if args.strandedness != "auto" else "unstranded"

        if layout == "PAIRED" and len(local_files) >= 2:
            # Sort to get R1/R2 order (_1 before _2)
            local_files.sort()
            rows.append({"id": safe_id, "fastq_1": local_files[0], "fastq_2": local_files[1], "strandedness": strand})
        else:
            rows.append({"id": safe_id, "fastq_1": local_files[0], "fastq_2": "", "strandedness": strand})

    if not rows:
        print("ERROR: no reads could be downloaded", file=sys.stderr)
        sys.exit(1)

    # Write sample sheet
    sample_sheet = outdir / "sample_sheet.csv"
    with open(sample_sheet, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Sample sheet written: {sample_sheet} ({len(rows)} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
