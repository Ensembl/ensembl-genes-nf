from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("GBS1", "unused")
os.environ.setdefault("GBP1", "3306")

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "bin"),
)

import resolve_species_inputs
from stable_id_mapping.ids import make_allocator, parse_id_range

import pytest


class FakeCursor:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def execute(self, query, parameters=None):
        pass

    def fetchone(self):
        return {
            "prefix": "ENSX",
            "range_start": 123,
            "range_end": 456,
        }


class FakeConnection:
    def cursor(self):
        return FakeCursor()

    def close(self):
        pass


def test_registry_ranges_are_formatted_to_exactly_11_digits(monkeypatch):
    monkeypatch.setattr(
        resolve_species_inputs,
        "_connect_ro",
        lambda: FakeConnection(),
    )

    ranges = resolve_species_inputs.resolve_stable_id_ranges("GCA_000000001.1")

    assert ranges == {
        "gene": "ENSXG:00000000123-00000000456",
        "transcript": "ENSXT:00000000123-00000000456",
        "translation": "ENSXP:00000000123-00000000456",
    }

    gene_range = parse_id_range(ranges["gene"])
    assert gene_range.width == 11
    assert make_allocator(gene_range, set()).allocate() == "ENSXG00000000123"

def mock_common_resolution(monkeypatch) -> None:
    monkeypatch.setattr(
        resolve_species_inputs,
        "fetch_core_db_meta",
        lambda db_name: {
            "assembly.accession": "GCA_000000001.3",
            "species.scientific_name": "Test species",
        },
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_target_assembly",
        lambda chain, version: (30, 300),
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_stable_id_ranges",
        lambda accession: {
            "gene": "ENSXG:00000000100-00000000199",
            "transcript": "ENSXT:00000000100-00000000199",
            "translation": "ENSXP:00000000100-00000000199",
        },
    )


def mock_nonconforming_stable_ids(monkeypatch) -> None:
    monkeypatch.setattr(
        resolve_species_inputs,
        "load_current_features",
        lambda db_name: {
            "gene": [{"stable_id": "WRONGG00000000100"}],
            "transcript": [{"stable_id": "WRONGT00000000100"}],
            "translation": [{"stable_id": "WRONGP00000000100"}],
            "exon": [{"stable_id": "WRONGE00000000100"}],
        },
    )


def test_reassign_does_not_resolve_reference_or_create_session(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)
    mock_nonconforming_stable_ids(monkeypatch)

    def unexpected_call(*args, **kwargs):
        raise AssertionError("mapping-only function was called")

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        unexpected_call,
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "insert_mapping_session",
        unexpected_call,
    )

    result = resolve_species_inputs.resolve_species_inputs(
        "test_core",
        requested_mode="reassign",
    )

    assert result.requested_mode == "reassign"
    assert result.effective_mode == "reassign"
    assert result.reference_assembly_id is None
    assert result.ref_fasta is None
    assert result.ref_gff is None
    assert result.mapping_session_id is None
    assert result.target_fasta is None
    assert result.target_gff is None


def test_auto_without_live_reference_uses_reassign(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)
    mock_nonconforming_stable_ids(monkeypatch)

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        lambda chain, version: None,
    )

    def unexpected_insert(*args, **kwargs):
        raise AssertionError("mapping session was created")

    monkeypatch.setattr(
        resolve_species_inputs,
        "insert_mapping_session",
        unexpected_insert,
    )

    result = resolve_species_inputs.resolve_species_inputs(
        "test_core",
        requested_mode="auto",
        target_fasta_override=Path("target.fa"),
        target_gff_override=Path("target.gff3"),
    )

    assert result.effective_mode == "reassign"
    assert result.mapping_session_id is None
    assert result.target_fasta is None
    assert result.target_gff is None


def test_auto_without_live_reference_and_conforming_ids_uses_no_action(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        lambda chain, version: None,
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "load_current_features",
        lambda db_name: {
            "gene": [{"stable_id": "ENSXG00000000100"}],
            "transcript": [{"stable_id": "ENSXT00000000100"}],
            "translation": [{"stable_id": "ENSXP00000000100"}],
            "exon": [{"stable_id": "ENSXE00000000100"}],
        },
    )

    def unexpected_insert(*args, **kwargs):
        raise AssertionError("mapping session was created")

    monkeypatch.setattr(
        resolve_species_inputs,
        "insert_mapping_session",
        unexpected_insert,
    )

    result = resolve_species_inputs.resolve_species_inputs(
        "test_core",
        requested_mode="auto",
    )

    assert result.effective_mode == "no_action"
    assert result.reference_assembly_id is None
    assert result.mapping_session_id is None
    assert result.ref_fasta is None
    assert result.ref_gff is None
    assert result.target_fasta is None
    assert result.target_gff is None


def test_auto_with_live_reference_uses_map_and_creates_session(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        lambda chain, version: (
            "GCA_000000001",
            2,
            20,
        ),
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "insert_mapping_session",
        lambda genebuild_status_id, reference_assembly_id: 777,
    )

    result = resolve_species_inputs.resolve_species_inputs(
        "test_core",
        requested_mode="auto",
        target_fasta_override=Path("target.fa"),
        target_gff_override=Path("target.gff3"),
        ref_fasta_override=Path("reference.fa"),
        ref_gff_override=Path("reference.gff3"),
    )

    assert result.effective_mode == "map"
    assert result.reference_assembly_id == 20
    assert result.mapping_session_id == 777
    assert result.ref_fasta == "reference.fa"
    assert result.ref_gff == "reference.gff3"


def test_explicit_map_without_live_reference_fails(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        lambda chain, version: None,
    )

    with pytest.raises(
        resolve_species_inputs.NoLiveReferenceError,
        match="Mode 'map' requested",
    ):
        resolve_species_inputs.resolve_species_inputs(
            "test_core",
            requested_mode="map",
            target_fasta_override=Path("target.fa"),
            target_gff_override=Path("target.gff3"),
        )

def test_resolve_pre_release_paths(
    monkeypatch,
    tmp_path: Path,
) -> None:
    pre_release_root = tmp_path / "pre-release"
    assembly_dir = (
        pre_release_root
        / "Larinus_planus"
        / "GCA_977012605.1"
    )
    assembly_dir.mkdir(parents=True)

    expected_fasta = (
        assembly_dir
        / "larinus_planus_gca977012605v1.dna.softmasked.fa.gz"
    )
    expected_gff = (
        assembly_dir
        / "larinus_planus_gca977012605v1.gff3.gz"
    )
    expected_fasta.touch()
    expected_gff.touch()

    monkeypatch.setattr(
        resolve_species_inputs,
        "PRE_RELEASE_ROOT",
        pre_release_root,
    )

    fasta, gff = resolve_species_inputs.resolve_pre_release_paths(
        "Larinus planus",
        "GCA_977012605",
        1,
    )

    assert fasta == expected_fasta
    assert gff == expected_gff


def test_mapping_falls_back_to_pre_release_for_target_only(
    monkeypatch,
) -> None:
    mock_common_resolution(monkeypatch)

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_live_reference",
        lambda chain, version: (
            "GCA_000000001",
            2,
            20,
        ),
    )

    def resolve_standard_paths(binomial_name, gca_chain, gca_version):
        if gca_version == 3:
            raise FileNotFoundError("standard target files missing")

        return (
            Path("reference.fa"),
            Path("reference.gff3"),
        )

    pre_release_calls = []

    def resolve_pre_release_paths(
        binomial_name,
        gca_chain,
        gca_version,
    ):
        pre_release_calls.append(
            (binomial_name, gca_chain, gca_version)
        )
        return (
            Path("pre-release-target.fa.gz"),
            Path("pre-release-target.gff3.gz"),
        )

    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_ftp_paths",
        resolve_standard_paths,
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "resolve_pre_release_paths",
        resolve_pre_release_paths,
    )
    monkeypatch.setattr(
        resolve_species_inputs,
        "insert_mapping_session",
        lambda genebuild_status_id, reference_assembly_id: 777,
    )

    result = resolve_species_inputs.resolve_species_inputs(
        "test_core",
        requested_mode="auto",
    )

    assert result.target_fasta == "pre-release-target.fa.gz"
    assert result.target_gff == "pre-release-target.gff3.gz"
    assert result.ref_fasta == "reference.fa"
    assert result.ref_gff == "reference.gff3"
    assert pre_release_calls == [
        ("Test species", "GCA_000000001", 3)
    ]
