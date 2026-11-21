import csv
import gzip
import hashlib
import io
import json
import sqlite3
import subprocess
import sys
from contextlib import closing
from pathlib import Path

import pytest


RIBOSEQ = Path(__file__).resolve().parents[1]
ACCESSION = "SRR11798519"
PROJECT = "PRJNA633140"


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend(str(RIBOSEQ / "bin"))
    import fetch_ena_run
    import prepare_rdp_processing

    return fetch_ena_run, prepare_rdp_processing


def table(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def test_rdp_export_uses_curated_type_and_keeps_held_runs(tmp_path, modules):
    _, exporter = modules
    db = tmp_path / "rdp snapshot.sqlite"
    rows = [
        (ACCESSION, "10090", "Mus musculus", "Ribo-Seq", "SINGLE"),
        ("SRR11798520", "10090", "Mus musculus", "Ribo-Seq", "PAIRED"),
        ("SRR11798521", "10090", "Mus musculus", "RNA-Seq", "SINGLE"),
        ("SRR11798522", "10090", "Mus musculus domesticus", "Ribo-Seq", "SINGLE"),
        ("SRR11798523", "10090", "Mus musculus", "Ribo-Seq", "SINGLE"),
    ]
    with closing(sqlite3.connect(db)) as conn:
        conn.execute(f"CREATE TABLE main_sample ({', '.join(exporter.SAMPLE_COLUMNS)})")
        conn.executemany(
            "INSERT INTO main_sample (Run, TaxID, ScientificName, LIBRARYTYPE, LibraryLayout, "
            "BioProject_id, LibraryStrategy) VALUES (?, ?, ?, ?, ?, ?, 'OTHER')",
            [(*row, PROJECT) for row in rows],
        )
        conn.commit()
    plan = tmp_path / "pilot.json"
    plan.write_text(json.dumps({
        "seed": 8, "quotas": {"mus_musculus": 1}, "controls": {ACCESSION: "exact run"}
    }))
    review = tmp_path / "review.csv"
    review.write_text("Run,reasons\nSRR11798523,study-level inference\n")
    out = tmp_path / "bundle"
    command = [
        sys.executable, str(RIBOSEQ / "bin/prepare_rdp_processing.py"), "--db", str(db),
        "--pilot-plan", str(plan), "--review-queue", str(review), "--out-dir", str(out),
    ]
    subprocess.run(command, check=True)
    sheet = [{"Run": ACCESSION, "study_accession": PROJECT}]
    example = RIBOSEQ / "species_run_setup/sample_sheets/mus_musculus_samplesheet.csv"
    header = example.read_bytes().splitlines(keepends=True)[0]
    for name in ("sample_sheets", "pilot", "controls"):
        path = out / name / example.name
        assert path.read_bytes() == header + f"{ACCESSION},{PROJECT}\r\n".encode()
        assert table(path, ",") == sheet
    inventory = {row["Run"]: row for row in table(out / "run_inventory.tsv")}
    assert len(inventory) == len(rows)
    assert inventory[ACCESSION]["raw_input_status"] == "unresolved"
    assert inventory["SRR11798520"]["hold_reasons"] == "outside_single_end_baseline"
    assert inventory["SRR11798521"]["hold_reasons"] == "outside_baseline_library_type"
    assert inventory["SRR11798522"]["reference_status"] == "unmapped_taxon"
    assert inventory["SRR11798523"]["hold_reasons"] == "rdp_metadata_review"
    reference = json.loads((out / "organism_setup/mus_musculus.json").read_text())
    assert reference["genome_gtf_url"].endswith("/GCA/000/001/635/9/ensembl/2026_04/geneset/genes.gtf.gz")
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 1
    assert "Output already exists" in result.stderr


def report_row(run=ACCESSION, layout="SINGLE", file_run=None):
    file_run = file_run or run
    return f"{run}\t{PROJECT}\tSRX8349996\t{layout}\t" \
           f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/019/{file_run}/{file_run}.fastq.gz\t{'a' * 32}"


def test_ena_selects_exact_run_even_if_the_report_has_siblings(modules):
    fetch, _ = modules
    report = "\n".join(("\t".join(fetch.FIELDS), report_row("SRR11798518"), report_row(),
                        report_row("SRR11798520"))) + "\n"
    result = fetch.parse_source(report, ACCESSION)
    assert result["run_accession"] == ACCESSION
    assert result["fastq_url"].endswith(f"/{ACCESSION}/{ACCESSION}.fastq.gz")


@pytest.mark.parametrize(("row", "reason"), [
    (report_row(run="SRR11798518"), "Found 0"),
    (report_row(layout="PAIRED"), "Unsupported ENA layout"),
    (report_row(file_run="SRR11798518"), "path does not identify the exact run"),
])
def test_ena_rejects_wrong_run_or_layout(modules, row, reason):
    fetch, _ = modules
    with pytest.raises(ValueError, match=reason):
        fetch.parse_source("\t".join(fetch.FIELDS) + "\n" + row, ACCESSION)


def test_ena_download_checks_archive_and_local_digests(modules, monkeypatch, tmp_path):
    fetch, _ = modules
    data = gzip.compress(b"@r1\nACGT\n+\nIIII\n")
    monkeypatch.setattr(fetch, "urlopen", lambda url, timeout: io.BytesIO(data))
    record = {"fastq_url": "https://ftp.sra.ebi.ac.uk/run.fastq.gz",
              "archive_md5": hashlib.md5(data, usedforsecurity=False).hexdigest()}
    path = tmp_path / "test.fastq.gz"
    assert fetch.download(record, path) == hashlib.sha256(data).hexdigest()
    assert path.read_bytes() == data
    with pytest.raises(ValueError, match="MD5 mismatch"):
        fetch.download({**record, "archive_md5": "0" * 32}, path)


def test_pre_release_reference_urls_are_preserved(modules):
    _, exporter = modules
    ref = {
        "beta_assembly": "GCA_052040795.2",
        "beta_provider": "pre_release",
        "beta_geneset_date": "2026-10-06",
        "genome_fasta_url": "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/Danio_rerio/GCA_052040795.2/danio_rerio_gca052040795v2.dna.softmasked.fa.gz",
        "gtf_url": "https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/Danio_rerio/GCA_052040795.2/danio_rerio_gca052040795v2.gtf.gz",
    }
    exporter.migrate_ensembl_urls(ref)
    assert ref["genome_fasta_url"].endswith("danio_rerio_gca052040795v2.dna.softmasked.fa.gz")
    assert ref["gtf_url"].endswith("danio_rerio_gca052040795v2.gtf.gz")
