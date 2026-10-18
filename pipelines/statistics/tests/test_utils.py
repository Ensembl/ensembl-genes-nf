import importlib.util
import io
import json
import zipfile
from pathlib import Path


def load_bin_script(script_name):
    script_path = Path(__file__).resolve().parents[1] / "bin" / script_name
    module_name = script_name[:-3]
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ena_assembly_path_builds_grouped_gca_url():
    fetch_genome = load_bin_script("fetch_genome.py")

    ena_base = "https://ftp.ebi.ac.uk/pub/databases/ena/assembly"

    assert (
        fetch_genome.ena_assembly_path(ena_base, "GCA_123456789.1")
        == "https://ftp.ebi.ac.uk/pub/databases/ena/assembly/GCA/123/456/789/GCA_123456789.1"
    )


def test_noncurrent_package_returns_current_accession():
    fetch_genome = load_bin_script("fetch_genome.py")
    for status in ("suppressed", "previous"):
        report = json.dumps(
            {
                "assemblyInfo": {"assemblyStatus": status},
                "currentAccession": "GCA_123456789.2",
            }
        ).encode()

        with io.BytesIO() as archive_buffer:
            with zipfile.ZipFile(archive_buffer, "w") as archive:
                archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl", report)
            with zipfile.ZipFile(io.BytesIO(archive_buffer.getvalue())) as archive:
                assert (
                    fetch_genome.suppressed_current_accession(archive, "GCA_123456789.1")
                    == "GCA_123456789.2"
                )
