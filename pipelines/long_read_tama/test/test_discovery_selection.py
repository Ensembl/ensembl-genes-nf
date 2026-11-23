from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "bin" / "discover_long_read_data.py"
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location("discover_long_read_data", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(run_id: str, sample: str, representation: str, tissue: str) -> dict:
    return {
        "run_accession": run_id,
        "sample_accession": sample,
        "representation_class": representation,
        "normalised_tissue": tissue,
        "developmental_stage": "adult",
        "base_count": "1000",
        "read_count": "100",
        "fastq_bytes": "1000",
    }


def test_same_sample_prefers_pacbio_subreads_over_ont():
    rows = [
        row("ONT1", "S1", "ONT_FASTQ", "brain"),
        row("PB1", "S1", "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE", "brain"),
    ]

    _inventory, selected = MODULE.select(rows, target=1, soft_download_budget=250.0, soft_raw_subread_budget=0)

    assert [item["run_accession"] for item in selected] == ["PB1"]
    assert selected[0]["platform_priority"] == 2
    assert selected[0]["selection_status"] == "SELECTED_WITH_WARNING"


def test_platform_tier_precedes_diversity_score():
    rows = [
        row("ONT1", "S1", "ONT_FASTQ", "brain"),
        row("PB1", "S2", "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE", "brain"),
    ]

    _inventory, selected = MODULE.select(rows, target=2, soft_download_budget=250.0, soft_raw_subread_budget=0)

    assert [item["run_accession"] for item in selected] == ["PB1", "ONT1"]
    assert selected[0]["selection_reason"].startswith("Platform tier 2")


def test_submitted_hifi_filename_is_consensus_evidence():
    assert MODULE._submitted_hifi_evidence(
        "ftp.sra.ebi.ac.uk/vol1/run/ERR1/movie.hifi_reads.fastq.gz"
    )
    assert MODULE._submitted_hifi_evidence(
        "ftp.sra.ebi.ac.uk/vol1/run/ERR1/movie.ccs.bam"
    )
    assert not MODULE._submitted_hifi_evidence(
        "ftp.sra.ebi.ac.uk/vol1/run/ERR1/movie.subreads.bam"
    )
