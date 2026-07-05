#!/usr/bin/env python3
"""Fetch genome fasta files from NCBI or ENA based on GCA accession."""
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import argparse
import os
import gzip
from pathlib import Path
import sys
import urllib.request
import zipfile
import shutil
import re
from urllib.error import HTTPError
import requests


def download_ncbi_assembly_report(gca: str, dest_folder: str| Path = ".") -> Path:
    """Download the assembly report from NCBI for a given GCA accession."""
    num = gca.split("_")[1].split(".")[0]

    p1, p2, p3 = num[0:3], num[3:6], num[6:9]

    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{p1}/{p2}/{p3}"

    r = requests.get(base, timeout=30)
    r.raise_for_status()

    match = re.search(rf"{gca}_[^\"/]+", r.text)

    if match is None:
        raise RuntimeError(f"Assembly directory not found for {gca}")

    assembly_dir = match.group(0)

    report_url = f"{base}/{assembly_dir}/{assembly_dir}_assembly_report.txt"

    dest_folder = Path(dest_folder)
    dest_folder.mkdir(parents=True, exist_ok=True)

    out_path = dest_folder / f"{assembly_dir}_assembly_report.txt"

    if out_path.exists():
        return out_path

    r = requests.get(report_url, timeout=30)
    r.raise_for_status()

    out_path.write_bytes(r.content)

    return out_path


def download_and_extract(url: str, output_dir: str| Path = ".") -> bool:
    """Download genome zip from NCBI and extract .fna files to output_dir.
    Returns True if successful, False otherwise.
    Inputs:
        url: URL to download the genome zip.
        output_dir: Directory to extract .fna files to.
    Returns:
        True if download and extraction were successful, False otherwise.
    """
    output_dir = Path(output_dir)
    zip_path = output_dir / "genome.zip"

    req = urllib.request.Request(url, headers={"Accept": "application/zip"})
    with urllib.request.urlopen(req) as r, open(zip_path, "wb") as out:
        out.write(r.read())
    extracted_roots: set[str] = set()
    with zipfile.ZipFile(zip_path) as z:
        files_to_extract = [f for f in z.namelist() if f.endswith(".fna")]
        if not any(f.endswith(".fna") for f in files_to_extract):
            return False
        for f in files_to_extract:
            z.extract(f, output_dir)

            src = output_dir / f
            dst = output_dir / Path(f).name

            shutil.move(src, dst)

            extracted_roots.add(Path(f).parts[0])
    # Cleanup extracted directory trees
    for root in extracted_roots:
        path = output_dir / root
        if path.is_dir():
            shutil.rmtree(path)

    zip_path.unlink()

    return True


def ena_assembly_path(ena_base: str, gca: str) -> str:
    """Construct ENA assembly path from base URL and GCA accession.
    Inputs:
        ena_base: Base URL for ENA assembly FTP.
        gca: GCA accession string.
    Returns:
        Full URL path to the assembly directory."""
    acc = gca.replace("GCA_", "")
    return f"{ena_base}/GCA/" f"{acc[0:3]}/{acc[3:6]}/{acc[6:9]}/{gca}"


def download_from_ena(ena_base: str, gca: str, output_dir: str) -> bool:
    """Download genome fasta from ENA for given GCA accession.
    Returns True if successful, False otherwise.
    Inputs:
        ena_base: Base URL for ENA assembly FTP.
        gca: GCA accession string.
        output_dir: Directory to save the downloaded fasta.
    Returns:
        True if download and extraction were successful, False otherwise.
    """
    try:
        base_url = ena_assembly_path(ena_base, gca)
        listing_url = f"{base_url}/"

        with urllib.request.urlopen(listing_url) as r:
            html = r.read().decode()

        # Find genomic fasta
        matches = [
            line.split('"')[1]
            for line in html.splitlines()
            if ".fna.gz" in line and "genomic" in line
        ]

        if not matches:
            return False

        fasta_gz = matches[0]
        fasta_url = f"{base_url}/{fasta_gz}"
        gz_path = os.path.join(output_dir, fasta_gz)

        print(f"Downloading genome from ENA: {fasta_url}")

        urllib.request.urlretrieve(fasta_url, gz_path)

        # Unzip
        fasta_path = gz_path[:-3]
        with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_path)
        return True

    except HTTPError:
        return False


def parse_assembly_report(assembly_report: Path) -> dict:
    """Build mapping GenBank accession -> Ensembl sequence name."""

    mapping = {}

    with open(assembly_report, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            seq_role = cols[1]
            assigned_molecule = cols[2]
            genbank = cols[4]

            if seq_role == "assembled-molecule":
                if assigned_molecule.lower() == "na":
                    mapping[genbank] = genbank
                else:
                    mapping[genbank] = assigned_molecule
            else:
                mapping[genbank] = genbank
    print(mapping)
    return mapping


def rewrite_fasta_headers(input_fasta: Path, output_fasta: Path, mapping: dict) -> None:
    """Rewrite FASTA headers using the assembly report mapping."""
    with open(input_fasta, "r", encoding="utf-8") as fin, open(
        output_fasta, "w", encoding="utf-8"
    ) as fout:
        for line in fin:
            if line.startswith(">"):
                accession = line[1:].split()[0]
                print(accession)
                fout.write(f">{mapping.get(accession, accession)}\n")

            else:
                fout.write(line)


def main():
    """Main function to parse arguments and download genome."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="fetch_genome.py 1.0.0")
    parser.add_argument("--gca", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--ncbi_base",
        default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession",
    )
    parser.add_argument(
        "--ena_base", default="https://ftp.ebi.ac.uk/pub/databases/ena/assembly"
    )
    parser.add_argument(
        "--reheader_file",
        action="store_true",
        help="Rewrite FASTA headers using the assembly report.",
    )
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    ncbi_url = (
        f"{args.ncbi_base}/{args.gca}/download"
        "?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
    )

    print(f"Downloading genome for {args.gca} from NCBI")

    if download_and_extract(ncbi_url, args.output_dir):
        print("Genome downloaded from NCBI")
        if args.reheader_file:
            output_dir = Path(args.output_dir)
            assembly_report = download_ncbi_assembly_report(args.gca, args.output_dir)
            input_fasta = next(output_dir.glob("*.fna"))
            output_fasta = output_dir / f"{input_fasta.stem}.fa"
            mapping = parse_assembly_report(assembly_report)
            rewrite_fasta_headers(input_fasta, output_fasta, mapping)
        sys.exit(0)
    print("NCBI genome not found, trying ENA")

    if download_from_ena(args.ena_base, args.gca, args.output_dir):
        print("Genome downloaded from ENA")
        sys.exit(0)

    print("Genome not found in NCBI or ENA", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
