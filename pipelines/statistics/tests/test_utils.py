import importlib.util
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
