#!/usr/bin/env python3
"""Resolve and verify one single-end FASTQ by its ENA run accession."""

import argparse
import csv
import hashlib
import io
import re
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from urllib.parse import urlencode, urlsplit
from urllib.request import urlopen


FILE_REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_HOST = "ftp.sra.ebi.ac.uk"
FIELDS = [
    "run_accession", "study_accession", "experiment_accession",
    "library_layout", "fastq_ftp", "fastq_md5",
]
SOURCE_FIELDS = [
    "run_accession", "study_accession", "experiment_accession", "library_layout",
    "fastq_url", "archive_md5",
]


def parse_source(report, accession):
    reader = csv.DictReader(io.StringIO(report), delimiter="\t")
    if not set(FIELDS).issubset(reader.fieldnames or []):
        raise ValueError("ENA file report is missing required columns")
    matches = [row for row in reader if row["run_accession"] == accession]
    if len(matches) != 1:
        raise ValueError(f"Expected one exact ENA record. Found {len(matches)}")
    row = matches[0]
    if row["library_layout"] != "SINGLE":
        raise ValueError(f"Unsupported ENA layout {row['library_layout']!r}. Expected SINGLE")
    ftp = row["fastq_ftp"] or ""
    md5 = row["fastq_md5"] or ""
    if not ftp or not md5:
        raise ValueError("ENA has no FASTQ URL or archive MD5")
    if ";" in ftp or ";" in md5:
        raise ValueError("Expected one FASTQ URL and one archive MD5 for a single-end run")
    url = urlsplit("https://" + ftp)
    path = PurePosixPath(url.path)
    if (url.netloc != ENA_HOST or url.query or url.fragment
            or path.parent.name != accession
            or path.name not in (f"{accession}.fastq.gz", f"{accession}_1.fastq.gz")):
        raise ValueError(f"ENA FASTQ path does not identify the exact run: {ftp!r}")
    if not re.fullmatch(r"[a-fA-F0-9]{32}", md5):
        raise ValueError(f"Invalid ENA FASTQ MD5 {md5!r}")
    return {
        **{key: row[key] for key in SOURCE_FIELDS[:4]},
        "fastq_url": url.geturl(),
        "archive_md5": md5.lower(),
    }


def resolve_source(accession):
    query = urlencode({
        "accession": accession, "result": "read_run",
        "fields": ",".join(FIELDS), "format": "tsv",
    })
    with urlopen(f"{FILE_REPORT}?{query}", timeout=60) as response:
        return parse_source(response.read().decode("utf-8"), accession)


def download(source, target):
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with urlopen(source["fastq_url"], timeout=120) as response, target.open("wb") as output:
        for chunk in iter(lambda: response.read(1024 * 1024), b""):
            output.write(chunk)
            md5.update(chunk)
            sha256.update(chunk)
    if md5.hexdigest() != source["archive_md5"]:
        raise ValueError(f"FASTQ MD5 mismatch. Expected {source['archive_md5']}, got {md5.hexdigest()}")
    return sha256.hexdigest()


def write_source(handle, record):
    writer = csv.DictWriter(handle, list(record), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerow(record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accession", help="SRR, ERR or DRR run accession")
    parser.add_argument("--metadata-only", action="store_true", help="Write the resolved ENA record to stdout")
    args = parser.parse_args()
    accession = args.accession
    if not re.fullmatch(r"[SED]RR[0-9]+", accession):
        parser.error("Expected an SRR, ERR or DRR run accession")
    try:
        source = resolve_source(accession)
        if args.metadata_only:
            write_source(sys.stdout, source)
            return
        fastq = Path(f"{accession}.fastq.gz")
        audit = Path(f"{accession}.fastq_source.tsv")
        with TemporaryDirectory(prefix=".ena-", dir=".") as directory:
            staged_fastq = Path(directory) / fastq
            staged_audit = Path(directory) / audit
            source["local_file"] = fastq.name
            source["local_sha256"] = download(source, staged_fastq)
            with staged_audit.open("w", newline="", encoding="utf-8") as handle:
                write_source(handle, source)
            staged_fastq.replace(fastq)
            staged_audit.replace(audit)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{accession}: {exc}\n")


if __name__ == "__main__":
    main()
