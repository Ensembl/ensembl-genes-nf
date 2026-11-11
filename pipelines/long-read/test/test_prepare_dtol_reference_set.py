from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "bin" / "prepare_dtol_reference_set.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_dtol_reference_set", MODULE_PATH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_reference_mapping_anchors_shared_contigs_before_renaming_primary_contigs():
    module = load_module()
    bam_dictionary = {
        "OZ023076.1": 100,
        "OZ023077.1": 90,
        "CAXDWV010000001.1": 50,
    }
    reference_dictionary = {
        "1": 100,
        "2": 90,
        "CAXDWV010000001.1": 50,
        "reference_only_contig": 10,
    }

    result = module.ordered_reference_mapping(bam_dictionary, reference_dictionary)

    assert result == (
        {
            "CAXDWV010000001.1": "CAXDWV010000001.1",
            "1": "OZ023076.1",
            "2": "OZ023077.1",
        },
        "renamed headers using shared contigs and unique contig lengths",
    )


def test_reference_mapping_rejects_ambiguous_length_only_rename():
    module = load_module()
    bam_dictionary = {"BAM_contig": 100}
    reference_dictionary = {"reference_one": 100, "reference_two": 100}

    assert module.ordered_reference_mapping(bam_dictionary, reference_dictionary) is None


def test_reference_dictionary_with_same_contigs_in_different_order_is_compatible():
    module = load_module()

    assert module.compare_dictionaries(
        {"bam_one": 100, "bam_two": 90},
        {"bam_two": 90, "bam_one": 100},
    ) == (
        "compatible",
        "BAM and reference contigs have the same names and lengths but a different order",
    )


def test_normalize_reference_accepts_renamed_primary_contigs_and_reordered_shared_contigs(tmp_path: Path):
    module = load_module()
    reference = tmp_path / "reference.fa"
    reference.write_text(
        ">1\nAAAAAAAAAA\n"
        ">2\nCCCCCCCCC\n"
        ">shared_two\nGGGG\n"
        ">shared_one\nTTT\n"
        ">reference_only\nAA\n"
    )
    output = tmp_path / "normalized.fa"
    mapping_path = tmp_path / "reference_header_map.tsv"
    bam_dictionary = {
        "OZ023076.1": 10,
        "OZ023077.1": 9,
        "shared_one": 3,
        "shared_two": 4,
    }

    normalized, detail = module.normalize_reference_to_bam(
        reference, bam_dictionary, output, mapping_path
    )

    assert normalized is True
    assert detail == "renamed headers using shared contigs and unique contig lengths"
    assert ">OZ023076.1" in output.read_text()
    assert ">OZ023077.1" in output.read_text()
    assert ">reference_only" not in output.read_text()
