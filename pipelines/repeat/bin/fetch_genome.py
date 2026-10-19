
#!/usr/bin/env python3

"""
Download genome assemblies from Ensembl Organisms or NCBI.

Priority:
    1. Ensembl Organisms
    2. NCBI Datasets API

Ensembl:
    - Discover the latest YYYY_MM release.
    - Download unmasked.fa.bgz.
    - Decompress.
    - Keep only the first token of FASTA headers.

NCBI:
    - Download genomic FASTA.
    - Download assembly report.
    - Reheader to Ensembl-style identifiers:
        assembled-molecule -> assigned molecule
        other sequences    -> GenBank accession

Final output:
    GCA_XXXXXXXXXvN.fa

Requirements:
    pip install requests
"""

import argparse
import gzip
import logging
import re
import shutil
import sys
import time
import zipfile

from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote

import requests


ENSEMBL_BASE = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"

NCBI_BASE = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/"
    "genome/accession"
)

NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

HEADERS = {
    "User-Agent": "GenomeDownloader/2.2",
}

LOG = logging.getLogger("fetch_genome")


# ============================================================
# HTTP utilities
# ============================================================

def get_with_retries(
    url: str,
    stream: bool = False,
    timeout: int = 60,
    max_retries: int = 5,
) -> requests.Response:
    """GET request with retries for transient errors."""

    for attempt in range(max_retries):
        try:
            response = requests.get(
                url,
                headers=HEADERS,
                stream=stream,
                timeout=timeout,
            )

            if response.status_code in (429, 500, 502, 503, 504):
                if attempt < max_retries - 1:
                    retry_after = response.headers.get("Retry-After")

                    try:
                        wait = (
                            min(float(retry_after), 300)
                            if retry_after
                            else min(2 ** attempt, 60)
                        )
                    except ValueError:
                        wait = min(2 ** attempt, 60)

                    response.close()
                    LOG.warning(
                        "HTTP error for %s. Retrying in %.1fs",
                        url,
                        wait,
                    )
                    time.sleep(wait)
                    continue

            response.raise_for_status()
            return response

        except (requests.Timeout, requests.ConnectionError):
            if attempt == max_retries - 1:
                raise

            wait = min(2 ** attempt, 60)
            LOG.warning("Connection error. Retrying in %ss", wait)
            time.sleep(wait)

    raise RuntimeError(f"Request failed: {url}")


def download_file(url: str, destination: Path) -> Path:
    """Download a file with streaming and atomic rename."""

    destination.parent.mkdir(parents=True, exist_ok=True)

    temporary = destination.with_name(
        destination.name + ".part"
    )

    try:
        with get_with_retries(
            url,
            stream=True,
            timeout=120,
        ) as response:

            with open(temporary, "wb") as fout:
                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):
                    if chunk:
                        fout.write(chunk)

        if temporary.stat().st_size == 0:
            raise RuntimeError(f"Downloaded empty file: {url}")

        temporary.replace(destination)

    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return destination


# ============================================================
# Directory listing parser
# ============================================================

class DirectoryParser(HTMLParser):
    """Extract directory links from an HTTP directory listing."""

    def __init__(self):
        super().__init__()
        self.directories = []

    def handle_starttag(self, tag, attrs):
        if tag != "a":
            return

        href = dict(attrs).get("href", "")

        if href.endswith("/") and href not in ("../", "./"):
            name = unquote(href.rstrip("/").split("/")[-1])
            self.directories.append(name)


def list_directories(url: str) -> list[str]:
    """List directories from an HTTP directory listing."""

    with get_with_retries(url, timeout=30) as response:
        parser = DirectoryParser()
        parser.feed(response.text)

    return parser.directories


# ============================================================
# Ensembl Organisms
# ============================================================

def ensembl_assembly_path(gca: str) -> str:
    """Construct Ensembl assembly path ending in /ensembl/."""

    match = re.fullmatch(r"GCA_(\d{9})\.(\d+)", gca)

    if not match:
        raise ValueError(f"Invalid GCA accession: {gca}")

    number, version = match.groups()

    return (
        f"{ENSEMBL_BASE}/GCA/"
        f"{number[:3]}/"
        f"{number[3:6]}/"
        f"{number[6:9]}/"
        f"{version}/ensembl/"
    )


def find_ensembl_genome(gca: str) -> str | None:
    """Find latest Ensembl release containing unmasked.fa.bgz."""

    base_url = ensembl_assembly_path(gca)

    LOG.info("Checking Ensembl: %s", base_url)

    try:
        releases = list_directories(base_url)

    except requests.HTTPError as exc:
        if (
            exc.response is not None
            and exc.response.status_code == 404
        ):
            LOG.info("%s not available in Ensembl", gca)
            return None
        raise

    releases = sorted(
        [
            r for r in releases
            if re.fullmatch(r"\d{4}_\d{2}", r)
            and 1 <= int(r.split("_")[1]) <= 12
        ],
        reverse=True,
    )

    if not releases:
        LOG.info("No Ensembl releases found for %s", gca)
        return None

    LOG.info("Available releases: %s", releases)

    for release in releases:

        genome_url = (
            f"{base_url}{release}/genome/unmasked.fa.bgz"
        )

        try:
            with get_with_retries(
                genome_url,
                stream=True,
                timeout=30,
            ):
                LOG.info("Selected Ensembl release: %s", release)
                return genome_url

        except requests.HTTPError as exc:
            if (
                exc.response is not None
                and exc.response.status_code == 404
            ):
                continue
            raise

    LOG.info("No Ensembl genome found for %s", gca)
    return None


def decompress_bgzf(
    compressed: Path,
    output_fasta: Path,
) -> Path:
    """Decompress BGZF to a temporary FASTA."""

    temporary = output_fasta.with_name(
        output_fasta.name + ".part"
    )

    try:
        with gzip.open(compressed, "rb") as fin:
            with open(temporary, "wb") as fout:
                shutil.copyfileobj(fin, fout)

        if temporary.stat().st_size == 0:
            raise RuntimeError("Decompressed Ensembl FASTA is empty")

        temporary.replace(output_fasta)

    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return output_fasta


def download_from_ensembl(
    gca: str,
    output_dir: Path,
) -> Path | None:
    """Download and decompress an Ensembl genome."""

    genome_url = find_ensembl_genome(gca)

    if genome_url is None:
        return None

    compressed = output_dir / f"{gca}.fa.bgz"
    raw_fasta = output_dir / f"{gca}.ensembl.raw.fa"

    LOG.info("Downloading Ensembl genome: %s", genome_url)

    try:
        download_file(genome_url, compressed)
        decompress_bgzf(compressed, raw_fasta)

    finally:
        compressed.unlink(missing_ok=True)

    return raw_fasta


# ============================================================
# NCBI assembly report
# ============================================================

def ncbi_assembly_base(gca: str) -> str:
    """Construct NCBI FTP parent directory."""

    match = re.fullmatch(r"GCA_(\d{9})\.(\d+)", gca)

    if not match:
        raise ValueError(f"Invalid GCA accession: {gca}")

    number = match.group(1)

    return (
        f"{NCBI_FTP}/GCA/"
        f"{number[:3]}/"
        f"{number[3:6]}/"
        f"{number[6:9]}"
    )


def download_ncbi_assembly_report(
    gca: str,
    dest_folder: str | Path = ".",
) -> Path:
    """Download NCBI assembly report."""

    dest_folder = Path(dest_folder)
    dest_folder.mkdir(parents=True, exist_ok=True)

    base = ncbi_assembly_base(gca)

    with get_with_retries(base) as response:
        html = response.text

    pattern = rf'({re.escape(gca)}_[^"/<>\s?]+)'
    matches = re.findall(pattern, html)

    if not matches:
        raise RuntimeError(
            f"Assembly directory not found for {gca}"
        )

    assembly_dir = matches[0]
    filename = f"{assembly_dir}_assembly_report.txt"

    report_url = f"{base}/{assembly_dir}/{filename}"
    out_path = dest_folder / filename

    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    LOG.info("Downloading assembly report: %s", report_url)

    return download_file(report_url, out_path)


# ============================================================
# NCBI Datasets genome download
# ============================================================

def download_and_extract(
    url: str,
    output_dir: Path,
) -> Path:
    """Download NCBI Datasets ZIP and extract genomic FASTA."""

    zip_path = output_dir / "genome.zip"

    download_file(url, zip_path)

    try:
        with zipfile.ZipFile(zip_path) as archive:

            members = [
                member
                for member in archive.infolist()
                if not member.is_dir()
                and member.filename.endswith(".fna")
            ]

            if len(members) != 1:
                raise RuntimeError(
                    "Expected one genomic FASTA in NCBI ZIP, "
                    f"found {len(members)}"
                )

            member = members[0]

            raw_fasta = output_dir / f"{member.filename.split('/')[-1]}"

            temporary = raw_fasta.with_name(
                raw_fasta.name + ".part"
            )

            try:
                with archive.open(member) as fin:
                    with open(temporary, "wb") as fout:
                        shutil.copyfileobj(fin, fout)

                if temporary.stat().st_size == 0:
                    raise RuntimeError("Extracted NCBI FASTA is empty")

                temporary.replace(raw_fasta)

            except Exception:
                temporary.unlink(missing_ok=True)
                raise

    finally:
        zip_path.unlink(missing_ok=True)

    return raw_fasta


def download_from_ncbi(
    gca: str,
    output_dir: Path,
) -> Path:
    """Download genomic FASTA using NCBI Datasets API."""

    ncbi_url = (
        f"{NCBI_BASE}/{gca}/download"
        "?include_annotation_type=GENOME_FASTA"
        "&hydrated=FULLY_HYDRATED"
    )

    LOG.info("Downloading genome from NCBI: %s", gca)

    return download_and_extract(
        ncbi_url,
        output_dir,
    )


# ============================================================
# Assembly report mapping
# ============================================================

def parse_assembly_report(
    assembly_report: Path,
) -> dict[str, str]:
    """
    Map NCBI accessions to Ensembl-style names.

    assembled-molecule -> assigned molecule
    other sequences    -> GenBank accession
    """

    mapping = {}

    with open(assembly_report, encoding="utf-8") as f:

        for line in f:

            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 7:
                raise ValueError(
                    f"Malformed assembly report row: {line}"
                )

            seq_role = cols[1]
            assigned_molecule = cols[2]
            genbank = cols[4]
            refseq = cols[6]

            if (
                seq_role == "assembled-molecule"
                and assigned_molecule.lower() != "na"
            ):
                ensembl_name = assigned_molecule
            else:
                ensembl_name = genbank

            if ensembl_name.lower() == "na":
                raise ValueError(
                    "Cannot determine name for "
                    f"assembly report row: {line.strip()}"
                )

            if genbank.lower() != "na":
                mapping[genbank] = ensembl_name

            if refseq.lower() != "na":
                mapping[refseq] = ensembl_name

    LOG.info("Loaded %d accession mappings", len(mapping))

    return mapping


# ============================================================
# FASTA header processing
# ============================================================

def process_fasta_headers(
    input_fasta: Path,
    output_fasta: Path,
    source: str,
    mapping: dict[str, str] | None = None,
) -> None:
    """
    Process FASTA headers safely.

    Ensembl:
        Keep first token only.

    NCBI:
        Replace accession with mapped Ensembl-style name.

    Always writes to a temporary file first.
    """

    temporary = output_fasta.with_name(
        output_fasta.name + ".part"
    )

    count = 0
    seen = set()

    try:
        with open(input_fasta, encoding="utf-8") as fin:
            with open(temporary, "w", encoding="utf-8") as fout:

                for line in fin:

                    if line.startswith(">"):
                        accession = line[1:].split()[0]

                        if source == "ensembl":
                            new_name = accession

                        elif source == "ncbi":
                            if mapping is None:
                                raise ValueError(
                                    "NCBI mapping is required"
                                )

                            if accession not in mapping:
                                raise ValueError(
                                    f"No mapping for {accession}"
                                )

                            new_name = mapping[accession]

                        else:
                            raise ValueError(
                                f"Unknown genome source: {source}"
                            )

                        if new_name in seen:
                            raise ValueError(
                                f"Duplicate FASTA header: {new_name}"
                            )

                        seen.add(new_name)
                        count += 1

                        fout.write(f">{new_name}\n")

                    else:
                        fout.write(line)

        if count == 0:
            raise RuntimeError(
                f"No FASTA headers found in {input_fasta}"
            )

        if temporary.stat().st_size == 0:
            raise RuntimeError(
                "Processed FASTA is empty"
            )

        temporary.replace(output_fasta)

    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    LOG.info(
        "Processed %d FASTA headers from %s",
        count,
        source,
    )


# ============================================================
# Main genome retrieval
# ============================================================

def fetch_genome(
    gca: str,
    output_dir: str | Path,
) -> Path:
    """
    Fetch genome from Ensembl or NCBI.

    Always download NCBI assembly report.
    Always produce cleaned FASTA.

    Final:
        GCA_XXXXXXXXXvN.fa
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gca_filename = gca.replace(".", "v")

    final_fasta = output_dir / f"{gca_filename}.fa"

    # Reuse existing nonempty final genome.
    if final_fasta.exists() and final_fasta.stat().st_size > 0:
        LOG.info("Using existing genome: %s", final_fasta)
        return final_fasta

    # --------------------------------------------------------
    # Step 1: Download assembly report (always)
    # --------------------------------------------------------

    assembly_report = download_ncbi_assembly_report(
        gca,
        output_dir,
    )

    # --------------------------------------------------------
    # Step 2: Try Ensembl first
    # --------------------------------------------------------

    raw_fasta = download_from_ensembl(
        gca,
        output_dir,
    )

    if raw_fasta is not None:
        source = "ensembl"
        mapping = None

    else:
        # ----------------------------------------------------
        # Step 3: NCBI fallback
        # ----------------------------------------------------

        source = "ncbi"

        raw_fasta = download_from_ncbi(
            gca,
            output_dir,
        )

        mapping = parse_assembly_report(
            assembly_report
        )

    # --------------------------------------------------------
    # Step 4: Process headers
    # --------------------------------------------------------

    process_fasta_headers(
        input_fasta=raw_fasta,
        output_fasta=final_fasta,
        source=source,
        mapping=mapping,
    )

    # --------------------------------------------------------
    # Step 5: Remove raw FASTA after successful processing
    # --------------------------------------------------------

    if raw_fasta != final_fasta:
        raw_fasta.unlink(missing_ok=True)

    LOG.info(
        "Final genome: %s (%s)",
        final_fasta,
        source,
    )

    return final_fasta


# ============================================================
# CLI
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Download genome from Ensembl or NCBI "
            "and standardize FASTA headers."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version="fetch_genome.py 2.2.0",
    )

    parser.add_argument(
        "--gca",
        required=True,
        help="Genome assembly accession, e.g. GCA_027574665.1",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    try:
        fasta = fetch_genome(
            gca=args.gca,
            output_dir=args.output_dir,
        )

        print(f"\nGenome ready: {fasta}")

    except Exception as exc:
        LOG.error("Genome retrieval failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
