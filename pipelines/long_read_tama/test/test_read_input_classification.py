import gzip
import importlib.util
from pathlib import Path
import subprocess

SPEC = importlib.util.spec_from_file_location("classification", Path(__file__).parents[1] / "bin" / "read_input_classification.py")
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def write_fastq(tmp_path, headers):
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt") as handle:
        for header in headers:
            handle.write(f"@{header}\nACGT\n+\n!!!!\n")
    return path


def test_subread_and_ccs_signatures():
    assert mod.header_signature("movie/42/10_20") == ("PACBIO_SUBREAD", "movie/42")
    assert mod.header_signature("movie/42/ccs") == ("PACBIO_CCS", "movie/42")


def test_uuid_requires_ont_evidence():
    uuid = "123e4567-e89b-12d3-a456-426614174000"
    assert mod.header_signature(uuid)[0] == "UNKNOWN"
    assert mod.header_signature(uuid, "@" + uuid + " runid=abc ch=4")[0] == "ONT"
    result = mod.classify(declared_platform="unknown", ena_platform="unknown", sra_platform="unavailable", artifacts=[{"artifact_role": "FASTQ"}], probe={"header_representation": "ONT", "ont_uuid_count": 1})
    assert result["classification"] == "UNKNOWN"


def test_probe_mixed_headers(tmp_path):
    result = mod.probe_fastq(write_fastq(tmp_path, ["movie/1/1_2", "movie/1/ccs"]))
    assert result["header_representation"] == "MIXED"
    assert result["records_sampled"] == 2


def test_subread_bam_conflict_is_quarantined():
    artifacts = mod.expand_artifacts("SRR1", "ENA", "https://example.org/a.subreads.bam", "d41d8cd98f00b204e9800998ecf8427e")
    result = mod.classify(declared_platform="PACBIO_SMRT", ena_platform="PACBIO_SMRT", sra_platform="ONT", artifacts=artifacts, probe={"header_representation": "ONT", "ont_uuid_count": 2})
    assert result["classification"] == "CONFLICT"
    assert result["status"] == "QUARANTINED"


def test_subread_fastq_is_quarantined():
    result = mod.classify(declared_platform="PACBIO_SMRT", ena_platform="PACBIO_SMRT", sra_platform="PACBIO", artifacts=[{"artifact_role": "FASTQ"}], probe={"header_representation": "PACBIO_SUBREAD", "subread_count": 2})
    assert result["classification"] == "PACBIO_SUBREAD_FASTQ_ONLY"
    assert result["proposed_action"] == "QUARANTINE"


def test_archive_subread_fastq_without_molecule_headers_is_quarantined():
    result = mod.classify(declared_platform="PACBIO_SMRT", ena_platform="PACBIO_SMRT", sra_platform="PACBIO",
                          library_strategy="RNA-Seq", library_source="TRANSCRIPTOMIC",
                          artifacts=[{"artifact_role": "FASTQ", "basename": "run_subreads.fastq.gz"}],
                          probe={"header_representation": "UNKNOWN"})
    assert result["classification"] == "PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE"
    assert result["confidence"] == "PROBABLE"


def test_ncbi_original_ccs_overrides_ena_subread_name():
    artifacts = [
        {"artifact_role": "FASTQ", "source": "ENA_FASTQ", "basename": "SRR26395028_subreads.fastq.gz"},
        {"artifact_role": "FASTQ", "source": "NCBI_ORIGINAL", "basename": "high.ccs.fq.gz"},
    ]
    result = mod.classify(declared_platform="PACBIO_SMRT", ena_platform="PACBIO_SMRT", sra_platform="PACBIO",
                          library_strategy="RNA-Seq", library_source="TRANSCRIPTOMIC",
                          artifacts=artifacts, probe={"header_representation": "UNKNOWN"})
    assert result["classification"] == "PACBIO_CCS_FASTQ"
    assert result["proposed_action"] == "ALIGN_PACBIO_CCS"
    assert result["confidence"] == "CONFIRMED"
    assert "NCBI_ORIGINAL_CCS_FASTQ" in result["reason_codes"]


def test_original_ccs_fastq_accepts_archive_rewritten_headers(tmp_path):
    path = write_fastq(tmp_path, ["SRR26395028.1", "SRR26395028.2"])
    result = mod.validate_fastq(path, "PACBIO_CCS_ORIGINAL")
    assert result["representation"] == "UNKNOWN"
    assert result["records"] == 2


def test_ncbi_original_processed_flnc_fastq_is_directly_usable():
    result = mod.classify(declared_platform="PACBIO_SMRT", ena_platform="PACBIO_SMRT", sra_platform="PACBIO",
                          artifacts=[{"artifact_role": "FASTQ", "source": "NCBI_ORIGINAL",
                                      "basename": "m84270_240911_210533_s4.skera.flnc.fastq.gz",
                                      "format": "PACBIO_PROCESSED_FASTQ"}],
                          probe={"header_representation": "UNKNOWN"})
    assert result["classification"] == "PACBIO_PROCESSED_FASTQ"
    assert result["proposed_action"] == "ALIGN_PACBIO_PROCESSED"


def test_artifact_lists_must_align():
    try:
        mod.expand_artifacts("SRR1", "ENA", "https://example.org/a.bam;https://example.org/a.bai", "d41d8cd98f00b204e9800998ecf8427e")
    except ValueError:
        return
    assert False, "mismatched artifact lists must fail"


def test_duplicate_molecules_fail_full_validation(tmp_path):
    path = write_fastq(tmp_path, ["movie/1/ccs", "movie/1/ccs"])
    try:
        mod.validate_fastq(path, "PACBIO_CCS")
    except ValueError as error:
        assert "duplicate" in str(error)
        return
    assert False, "duplicate CCS molecule IDs must fail"


def test_invalid_gzip_fails_full_validation(tmp_path):
    path = tmp_path / "corrupt.fastq.gz"
    path.write_bytes(b"not gzip")
    try:
        mod.validate_fastq(path, "ONT")
    except (OSError, EOFError):
        return
    assert False, "corrupt gzip must fail validation"


def test_incomplete_fastq_fails_full_validation(tmp_path):
    path = tmp_path / "incomplete.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@read-1\nACGT\n+\n")
    try:
        mod.validate_fastq(path, "ONT")
    except ValueError as error:
        assert "truncated" in str(error) or "complete" in str(error)
        return
    assert False, "incomplete FASTQ must fail validation"


def test_versioned_candidate_manifest_preserves_declared_platform(tmp_path):
    source = tmp_path / "candidate.tsv"
    source.write_text("run_accession\ttissue\tdescription\turl\tmd5\tplatform\nSRR1\tliver\tdesc\thttps://example.org/a.fastq.gz\td41d8cd98f00b204e9800998ecf8427e\tONT\n")
    output, report = tmp_path / "normalised.tsv", tmp_path / "report.tsv"
    subprocess.run(["python3", str(Path(__file__).parents[1] / "bin" / "normalise_manifest.py"), str(source), str(output), str(report)], check=True)
    assert "\tONT\t" in output.read_text()


def test_approved_gate_rejects_conflict():
    row = {"status": "APPROVED", "review_decision": "APPROVE", "classification": "CONFLICT", "proposed_action": "QUARANTINE", "run_accession": "SRR1", "selected_artifact_md5": "d41d8cd98f00b204e9800998ecf8427e", "selected_artifact_uri": "https://example.org/a.fastq.gz", "selected_artifact_basename": "a.fastq.gz", "minimap2_preset": "splice"}
    try:
        mod.validate_approved([row])
    except ValueError:
        return
    assert False, "conflict must not pass the approval gate"


def test_auto_approve_safe_filters_quarantined_rows(tmp_path):
    report_dir = tmp_path / "classification_report"
    report_dir.mkdir()
    (report_dir / "run_classification.tsv").write_text(
        "run_accession\ttissue\tdescription\tclassification\tproposed_action\tselected_artifact_uri\tselected_artifact_md5\tselected_artifact_basename\treason_codes\tstatus\n"
        "SRR1\tliver\tgood\tONT_FASTQ\tALIGN_ONT\thttps://example.org/a.fastq.gz\td41d8cd98f00b204e9800998ecf8427e\ta.fastq.gz\tOK\tREADY_FOR_REVIEW\n"
        "SRR2\tbrain\tbad\tCONFLICT\tQUARANTINE\thttps://example.org/b.fastq.gz\td41d8cd98f00b204e9800998ecf8427e\tb.fastq.gz\tPLATFORM_CONFLICT\tQUARANTINED\n"
    )
    output = tmp_path / "approved.tsv"
    audit = tmp_path / "audit.tsv"
    subprocess.run(["python3", str(Path(__file__).parents[1] / "bin" / "auto_approve_selected.py"), str(report_dir), str(output), str(audit)], check=True)
    assert "SRR1" in output.read_text()
    assert "SRR2" not in output.read_text()
    assert "SRR2\tCONFLICT\tQUARANTINED\tQUARANTINE" in audit.read_text()


def test_auto_approve_safe_reports_empty_selection(tmp_path):
    report_dir = tmp_path / "classification_report"
    report_dir.mkdir()
    (report_dir / "run_classification.tsv").write_text(
        "run_accession\tclassification\tstatus\treason_codes\n"
        "SRR2\tCONFLICT\tQUARANTINED\tPLATFORM_CONFLICT\n"
    )
    output = tmp_path / "approved.tsv"
    audit = tmp_path / "audit.tsv"
    result = subprocess.run(
        ["python3", str(Path(__file__).parents[1] / "bin" / "auto_approve_selected.py"),
         str(report_dir), str(output), str(audit)], capture_output=True, text=True
    )
    assert result.returncode != 0
    assert "no safely runnable runs" in result.stderr
    assert "SRR2\tCONFLICT\tQUARANTINED\tQUARANTINE" in audit.read_text()
