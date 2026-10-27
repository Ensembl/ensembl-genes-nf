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


BUSCO_GENOME_SUMMARY = """# BUSCO version is: 6.0.0
# The lineage dataset is: primates_odb12 (Creation date: 2024-01-08, number of genomes: 23, number of BUSCOs: 13780)
# Summarized benchmarking in BUSCO notation for file genome.fa
# BUSCO was run in mode: euk_genome_min

	***** Results: *****

	C:98.2%[S:96.1%,D:2.1%],F:0.9%,M:0.9%,n:13780,E:1.2%
	13532	Complete BUSCOs (C)
	13243	Complete and single-copy BUSCOs (S)
	289	Complete and duplicated BUSCOs (D)
	124	Fragmented BUSCOs (F)
	124	Missing BUSCOs (M)
	13780	Total BUSCO groups searched
"""


def test_busco_asmdb_parses_genome_summary(tmp_path):
    busco_asmdb = load_bin_script("busco_metakeys_to_asmdb.py")
    summary = tmp_path / "short_summary.txt"
    summary.write_text(BUSCO_GENOME_SUMMARY)

    metrics = busco_asmdb.parse_busco_file(str(summary))

    assert metrics["assembly.busco"] == "C:98.2%[S:96.1%,D:2.1%],F:0.9%,M:0.9%,n:13780,E:1.2%"
    assert metrics["assembly.busco_version"] == "6.0.0"
    assert metrics["assembly.busco_dataset"] == "primates_odb12"
    assert metrics["assembly.busco_mode"] == "genome"
    assert metrics["assembly.busco_completeness"] == "13532"
    assert metrics["assembly.busco_total"] == 13780
    assert metrics["assembly.busco_erroneus"] == "1.2"


def test_busco_asmdb_rejects_protein_summary(tmp_path):
    busco_asmdb = load_bin_script("busco_metakeys_to_asmdb.py")
    summary = tmp_path / "short_summary.txt"
    summary.write_text(BUSCO_GENOME_SUMMARY.replace("euk_genome_min", "proteins"))

    try:
        busco_asmdb.parse_busco_file(str(summary))
    except ValueError as err:
        assert "proteins" in str(err)
    else:
        raise AssertionError("Expected ValueError for protein mode summary")


def test_busco_asmdb_erroneous_from_score_line_regardless_of_mode(tmp_path):
    busco_asmdb = load_bin_script("busco_metakeys_to_asmdb.py")
    summary = tmp_path / "short_summary.txt"

    summary.write_text(BUSCO_GENOME_SUMMARY.replace("euk_genome_min", "genome"))
    assert busco_asmdb.parse_busco_file(str(summary))["assembly.busco_erroneus"] == "1.2"

    summary.write_text(BUSCO_GENOME_SUMMARY.replace(",E:1.2%", ""))
    metrics = busco_asmdb.parse_busco_file(str(summary))
    assert "assembly.busco_erroneus" not in metrics
    assert metrics["assembly.busco"] == "C:98.2%[S:96.1%,D:2.1%],F:0.9%,M:0.9%,n:13780"


def test_busco_asmdb_parses_connection_params_string():
    busco_asmdb = load_bin_script("busco_metakeys_to_asmdb.py")

    params = busco_asmdb.parse_connection_params(
        '{"host": "h", "port": "4527", "user": "u", "password": "p\'w", "database": "d"}'
    )
    assert params == {"host": "h", "port": 4527, "user": "u", "password": "p'w", "database": "d"}

    for bad in ("not json", '["a", "list"]'):
        try:
            busco_asmdb.parse_connection_params(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Expected ValueError for {bad!r}")
