#!/usr/bin/env python3
import argparse
import os
import gzip
import sys
import urllib.request
import zipfile
import shutil
from urllib.error import HTTPError

def download_and_extract(url: str, output_dir: str) -> bool:
    zip_path = os.path.join(output_dir, "genome.zip")

    req = urllib.request.Request(url, headers={"Accept": "application/zip"})
    with urllib.request.urlopen(req) as r, open(zip_path, "wb") as out:
        out.write(r.read())
    extracted_roots: set[str] = set()
    with zipfile.ZipFile(zip_path) as z:
        fna_files = [f for f in z.namelist() if f.endswith(".fna")]
        if not fna_files:
            return False

        for f in fna_files:
            z.extract(f, output_dir)
            shutil.move(
                os.path.join(output_dir, f),
                os.path.join(output_dir, os.path.basename(f))
            )
            # Track top-level extracted directory (e.g. ncbi_dataset)
            #extracted_roots.add(f.split(os.sep)[0])
            extracted_roots.add(f.split("/")[0])

        # Cleanup extracted directory trees
    for root in extracted_roots:
        path = os.path.join(output_dir, root)
        if os.path.isdir(path):
            shutil.rmtree(path)

    os.remove(zip_path)
    return True
def ena_assembly_path(ena_base: str, gca: str) -> str:
    acc = gca.replace("GCA_", "")
    return (
        f"{ena_base}/GCA/"
        f"{acc[0:3]}/{acc[3:6]}/{acc[6:9]}/{gca}"
    )
def download_from_ena(ena_base: str, gca: str, output_dir: str) -> bool:
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
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version",
        action="version",
        version="fetch_genome.py 1.0.0"
    )
    parser.add_argument("--gca", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--ncbi_base",
        default="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
    )
    parser.add_argument(
        "--ena_base",
        default="https://ftp.ebi.ac.uk/pub/databases/ena/assembly"
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
        sys.exit(0)
    print("NCBI genome not found, trying ENA")

    if download_from_ena(args.ena_base, args.gca, args.output_dir):
        print("Genome downloaded from ENA")
        sys.exit(0)

    print("Genome not found in NCBI or ENA", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()

