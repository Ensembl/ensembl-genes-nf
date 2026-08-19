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
import json
import sys
import urllib.request
import zipfile
import shutil
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError


def suppressed_current_accession(zip_file: zipfile.ZipFile, requested_accession: str) -> str | None:
    """Return NCBI's replacement accession for a non-current assembly package."""
    report_name = "ncbi_dataset/data/assembly_data_report.jsonl"
    try:
        report = json.loads(zip_file.read(report_name).decode().splitlines()[0])
    except (KeyError, IndexError, json.JSONDecodeError, UnicodeDecodeError):
        return None

    assembly_info = report.get("assemblyInfo", {})
    current_accession = report.get("currentAccession") or assembly_info.get("currentAccession")
    if (
        assembly_info.get("assemblyStatus") in {"suppressed", "previous"}
        and current_accession
        and current_accession != requested_accession
    ):
        return current_accession
    return None


def download_and_extract(
    url: str, output_dir: str, requested_accession: str, allow_suppressed: bool = False
) -> bool | str:
    """Download genome zip from NCBI and extract .fna files to output_dir.

    Returns True if successful, False otherwise, or a replacement accession
    when an explicitly enabled non-current accession has a currentAccession.
    Inputs:
        url: URL to download the genome zip.
        output_dir: Directory to extract .fna files to.
    Returns:
        True if download and extraction succeeded, False otherwise, or the
        current replacement accession for an allowed non-current assembly.
    """
    zip_path = os.path.join(output_dir, "genome.zip")

    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/zip"})
            with urllib.request.urlopen(req) as r, open(zip_path, "wb") as out:
                # Stream the archive so a large download does not occupy an
                # unnecessary second copy of the ZIP in Python memory.
                shutil.copyfileobj(r, out)

            extracted_roots: set[str] = set()
            with zipfile.ZipFile(zip_path) as z:
                fna_files = [f for f in z.namelist() if f.endswith(".fna")]
                if not fna_files:
                    if allow_suppressed:
                        replacement = suppressed_current_accession(z, requested_accession)
                        if replacement:
                            os.remove(zip_path)
                            return replacement
                    return False

                for f in fna_files:
                    z.extract(f, output_dir)
                    shutil.move(os.path.join(output_dir, f), os.path.join(output_dir, os.path.basename(f)))
                    extracted_roots.add(f.split("/")[0])

            for root in extracted_roots:
                path = os.path.join(output_dir, root)
                if os.path.isdir(path):
                    shutil.rmtree(path)

            os.remove(zip_path)
            return True
        except (HTTPError, URLError, IncompleteRead, OSError, zipfile.BadZipFile) as exc:
            if os.path.exists(zip_path):
                os.remove(zip_path)
            if attempt == 2:
                print(f"NCBI download failed after retries: {exc}", file=sys.stderr)
            else:
                time.sleep(15 * (attempt + 1))
    return False


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
            line.split('"')[1] for line in html.splitlines() if ".fna.gz" in line and "genomic" in line
        ]

        if not matches:
            return False

        fasta_gz = matches[0]
        fasta_url = f"{base_url}/{fasta_gz}"
        gz_path = os.path.join(output_dir, fasta_gz)

        print(f"Downloading genome from ENA: {fasta_url}")

        for attempt in range(3):
            try:
                with urllib.request.urlopen(fasta_url) as response, open(gz_path, "wb") as output:
                    shutil.copyfileobj(response, output)
                break
            except (HTTPError, URLError, IncompleteRead, OSError):
                if os.path.exists(gz_path):
                    os.remove(gz_path)
                if attempt == 2:
                    return False
                time.sleep(15 * (attempt + 1))

        # Unzip
        fasta_path = gz_path[:-3]
        with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_path)
        return True

    except (HTTPError, URLError, OSError):
        return False


def main():
    """Main function to parse arguments and download genome."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="fetch_genome.py 1.0.0")
    parser.add_argument("--gca", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--allow-suppressed-accessions",
        action="store_true",
        help="Use NCBI's currentAccession for suppressed or previous assemblies",
    )
    parser.add_argument(
        "--ncbi_base", default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
    )
    parser.add_argument("--ena_base", default="https://ftp.ebi.ac.uk/pub/databases/ena/assembly")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    ncbi_url = (
        f"{args.ncbi_base}/{args.gca}/download"
        "?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
    )

    print(f"Downloading genome for {args.gca} from NCBI")

    download_result = download_and_extract(
        ncbi_url, args.output_dir, args.gca, args.allow_suppressed_accessions
    )
    if download_result is True:
        print("Genome downloaded from NCBI")
        sys.exit(0)

    if isinstance(download_result, str):
        replacement_url = (
            f"{args.ncbi_base}/{download_result}/download"
            "?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED"
        )
        print(
            f"Requested assembly {args.gca} is non-current; "
            f"using current accession {download_result} for genome BUSCO",
            file=sys.stderr,
        )
        if download_and_extract(replacement_url, args.output_dir, download_result) is True:
            print(
                f"Genome downloaded from NCBI using {download_result} "
                f"(requested {args.gca})"
            )
            sys.exit(0)
    print("NCBI genome not found, trying ENA")

    if download_from_ena(args.ena_base, args.gca, args.output_dir):
        print("Genome downloaded from ENA")
        sys.exit(0)

    print("Genome not found in NCBI or ENA", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
