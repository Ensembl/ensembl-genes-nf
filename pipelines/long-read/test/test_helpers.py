from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[1]
BIN = ROOT / "bin"
DATA = ROOT / "test" / "data"


def run_helper(name: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(BIN / name), *args], text=True, capture_output=True)


def test_manifest_normalization_resolves_relative_input(tmp_path: Path):
    result = run_helper(
        "normalize_manifest.py",
        "--manifest", str(DATA / "fastq_manifest.tsv"),
        "--mode", "fastq",
        "--cohort-id", "fixture",
        "--normalized", str(tmp_path / "normalized.tsv"),
        "--manifest-checksum", str(tmp_path / "normalized.sha256"),
        "--checksums", str(tmp_path / "checksums.tsv"),
        "--status", str(tmp_path / "status.tsv"),
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader((tmp_path / "normalized.tsv").open(), delimiter="\t"))
    assert rows[0]["input"].endswith("/test/data/reads.fastq")
    assert len(rows[0]["input_sha256"]) == 64


def test_manifest_normalization_requires_and_records_bam_index(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    bam = tmp_path / "reads.bam"
    bai = tmp_path / "reads.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"index")
    manifest.write_text(
        "id\tcohort_id\trun_id\ttechnology\task\tassembly_id\tread_set_id\tsample_id\tlibrary_id\tread_representation\tsource_uri\tinput\n"
        f"run1\tfixture\trun1\tpacbio\tfixture\tfixture\tset1\tsample1\tlib1\tbam\t{bam.name}\t{bam.name}\n"
    )
    result = run_helper(
        "normalize_manifest.py", "--manifest", str(manifest), "--mode", "bam", "--cohort-id", "fixture",
        "--normalized", str(tmp_path / "normalized.tsv"), "--manifest-checksum", str(tmp_path / "normalized.sha256"),
        "--checksums", str(tmp_path / "checksums.tsv"), "--status", str(tmp_path / "status.tsv")
    )
    assert result.returncode == 0, result.stderr
    assert "[index]" in (tmp_path / "checksums.tsv").read_text()

    bai.unlink()
    failed = run_helper(
        "normalize_manifest.py", "--manifest", str(manifest), "--mode", "bam", "--cohort-id", "fixture",
        "--normalized", str(tmp_path / "normalized2.tsv"), "--manifest-checksum", str(tmp_path / "normalized2.sha256"),
        "--checksums", str(tmp_path / "checksums2.tsv"), "--status", str(tmp_path / "status2.tsv")
    )
    assert failed.returncode != 0


def test_manifest_normalization_requires_pbi_for_explicit_pbccs(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    bam = tmp_path / "subreads.bam"
    pbi = tmp_path / "subreads.bam.pbi"
    bam.write_bytes(b"bam")
    pbi.write_bytes(b"pbi")
    manifest.write_text(
        "id\tcohort_id\trun_id\ttechnology\tassembly_id\tread_set_id\tsample_id\tlibrary_id\tread_representation\tsource_uri\tinput\n"
        f"run1\tfixture\trun1\tpacbio\tfixture\tset1\tsample1\tlib1\traw_subreads\t{bam.name}\t{bam.name}\n"
    )
    result = run_helper(
        "normalize_manifest.py", "--manifest", str(manifest), "--mode", "bam", "--ccs-mode", "pbccs",
        "--cohort-id", "fixture", "--normalized", str(tmp_path / "normalized.tsv"),
        "--manifest-checksum", str(tmp_path / "normalized.sha256"), "--checksums", str(tmp_path / "checksums.tsv"),
        "--status", str(tmp_path / "status.tsv")
    )
    assert result.returncode == 0, result.stderr
    assert "subreads.bam.pbi" in (tmp_path / "normalized.tsv").read_text()


def test_retrospective_manifest_resolves_historical_bed12(tmp_path: Path):
    result = run_helper(
        "normalize_manifest.py", "--manifest", str(DATA / "retrospective_manifest.tsv"),
        "--mode", "retrospective", "--cohort-id", "fixture",
        "--normalized", str(tmp_path / "normalized.tsv"),
        "--manifest-checksum", str(tmp_path / "normalized.sha256"),
        "--checksums", str(tmp_path / "checksums.tsv"),
        "--status", str(tmp_path / "status.tsv"),
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader((tmp_path / "normalized.tsv").open(), delimiter="\t"))
    assert rows[0]["bed12"].endswith("/test/data/retrospective.bed12")
    assert len(rows[0]["bed12_sha256"]) == 64


def test_manifest_allows_remote_input_only_for_explicit_acquisition(tmp_path: Path):
    manifest = tmp_path / "remote.tsv"
    manifest.write_text(
        "id\tcohort_id\trun_id\ttechnology\task\tassembly_id\tread_set_id\tsample_id\tlibrary_id\tread_representation\tsource_uri\n"
        "run1\tfixture\trun1\tpacbio\tfixture\tfixture\tset1\tsample1\tlib1\tfastq\thttps://example.invalid/reads.fastq.gz\n"
    )
    result = run_helper(
        "normalize_manifest.py", "--manifest", str(manifest), "--mode", "approved_manifest",
        "--allow-unresolved-inputs", "--cohort-id", "fixture",
        "--normalized", str(tmp_path / "normalized.tsv"), "--manifest-checksum", str(tmp_path / "normalized.sha256"),
        "--checksums", str(tmp_path / "checksums.tsv"), "--status", str(tmp_path / "status.tsv")
    )
    assert result.returncode == 0, result.stderr
    row = next(csv.DictReader((tmp_path / "normalized.tsv").open(), delimiter="\t"))
    assert row["input"] == "https://example.invalid/reads.fastq.gz"
    assert row["input_sha256"] == ""


def test_finalize_prepared_manifest_replaces_input_checksum(tmp_path: Path):
    normalized = tmp_path / "normalized.tsv"
    normalized.write_text("id\tinput\tinput_sha256\nrun1\tremote\t\n")
    prepared = tmp_path / "prepared.tsv"
    prepared.write_text("id\tinput\tinput_sha256\tinput_index\tinput_index_sha256\nrun1\tstaged.fastq\t" + "a" * 64 + "\t\t\n")
    result = run_helper(
        "finalize_prepared_manifest.py", "--normalized", str(normalized), "--prepared", str(prepared),
        "--output", str(tmp_path / "final.tsv"), "--checksums", str(tmp_path / "checksums.tsv"),
        "--status", str(tmp_path / "status.tsv")
    )
    assert result.returncode == 0, result.stderr
    row = next(csv.DictReader((tmp_path / "final.tsv").open(), delimiter="\t"))
    assert row["input"] == "staged.fastq"
    assert row["input_sha256"] == "a" * 64


def test_validate_fastq_handles_gzip_without_materialising_input(tmp_path: Path):
    reads = tmp_path / "reads.fastq.gz"
    with gzip.open(reads, "wt") as handle:
        handle.write(DATA.joinpath("reads.fastq").read_text())
    status = tmp_path / "status.tsv"
    result = run_helper("validate_fastq.py", "--reads", str(reads), "--run-id", "run1", "--status", str(status))
    assert result.returncode == 0, result.stderr
    assert "run1\tVALID\t1" in status.read_text()


def test_merge_filelist_is_sorted_and_checksum_bearing(tmp_path: Path):
    first = tmp_path / "z.gtf"
    second = tmp_path / "a.gtf"
    first.write_text("z\n")
    second.write_text("a\n")
    file_list = tmp_path / "merge_filelist.txt"
    checksum = tmp_path / "merge_filelist.sha256"
    result = run_helper(
        "merge_filelist.py", "--inputs", str(first), str(second),
        "--file-list", str(file_list), "--checksum", str(checksum)
    )
    assert result.returncode == 0, result.stderr
    assert file_list.read_text().splitlines() == sorted([str(first.resolve()), str(second.resolve())])
    assert len(checksum.read_text().split()[0]) == 64


def test_merge_filelist_manifest_records_inputs_and_checksum(tmp_path: Path):
    first = tmp_path / "b.gtf"
    second = tmp_path / "a.gtf"
    first.write_text("b\n")
    second.write_text("a\n")
    file_list = tmp_path / "merge_filelist.txt"
    checksum = tmp_path / "merge_filelist.sha256"
    manifest = tmp_path / "merge_manifest.json"
    result = run_helper(
        "merge_filelist.py", "--inputs", str(first), str(second),
        "--file-list", str(file_list), "--checksum", str(checksum),
        "--manifest", str(manifest), "--backend", "tmerge",
        "--scope", "run", "--status", "COMPLETE"
    )
    assert result.returncode == 0, result.stderr
    record = json.loads(manifest.read_text())
    assert record["backend"] == "tmerge"
    assert record["scope"] == "run"
    assert record["status"] == "COMPLETE"
    assert record["inputs"] == sorted([str(first.resolve()), str(second.resolve())])
    assert record["file_list_sha256"] == checksum.read_text().split()[0]


def test_isoquant_manifest_records_content_checksums(tmp_path: Path):
    reference = tmp_path / "reference.fa"
    bam = tmp_path / "reads.bam"
    genedb = tmp_path / "genes.gtf"
    reference.write_text(">chr1\nAAAA\n")
    bam.write_bytes(b"bam")
    genedb.write_text("gene\n")
    result = run_helper(
        "write_isoquant_manifest.py", "--reference", str(reference), "--bam", str(bam),
        "--genedb", str(genedb), "--data-type", "pacbio_ccs", "--analysis", "transcript_discovery",
        "--large-output", "read_info", "read_assignments", "--args", "--threads 2",
        "--output", str(tmp_path / "isoquant.json")
    )
    assert result.returncode == 0, result.stderr
    record = json.loads((tmp_path / "isoquant.json").read_text())
    assert len(record["reference"]["sha256"]) == 64
    assert len(record["bam"]["sha256"]) == 64
    assert len(record["genedb"]["sha256"]) == 64
    assert record["large_output"] == ["read_info", "read_assignments"]


def test_flair_manifest_records_native_product_checksums(tmp_path: Path):
    reference = tmp_path / "reference.fa"
    native = tmp_path / "native_products"
    reference.write_text(">chr1\nAAAA\n")
    native.mkdir()
    (native / "flair.isoforms.gtf").write_text("gtf\n")
    (native / "flair.isoforms.bed").write_text("bed\n")
    result = run_helper(
        "write_flair_manifest.py", "--output", str(tmp_path / "flair.json"),
        "--native-products", str(native), "--run-id", "run1", "--shard-id", "run1__chr1",
        "--status", "SUCCESS", "--input", str(reference)
    )
    assert result.returncode == 0, result.stderr
    record = json.loads((tmp_path / "flair.json").read_text())
    assert record["run_id"] == "run1"
    assert {item["name"] for item in record["native_products"]} == {
        "flair.isoforms.bed", "flair.isoforms.gtf"
    }
    assert len(record["inputs"][0]["sha256"]) == 64


def test_container_smoke_discovers_all_pinned_tools(tmp_path: Path):
    result = run_helper(
        "container_smoke.py", str(ROOT / "resources" / "containers.yml"),
        "--engine", "/usr/bin/true", "--output", str(tmp_path / "smoke.json")
    )
    assert result.returncode == 0, result.stderr
    record = json.loads((tmp_path / "smoke.json").read_text())
    assert len(record["results"]) >= 11
    assert {item["name"] for item in record["results"]} >= {"tama", "isoquant", "flair", "bambu", "python_helper"}


def test_container_smoke_timeout_output_is_serializable(tmp_path: Path):
    helper = importlib.util.spec_from_file_location(
        "container_smoke", BIN / "container_smoke.py"
    )
    assert helper is not None and helper.loader is not None
    module = importlib.util.module_from_spec(helper)
    helper.loader.exec_module(module)
    assert module.text_value(b"timeout") == "timeout"
    assert module.text_value(None) == ""
    returncode, stdout, stderr = module.run_command(
        [sys.executable, "-c", "import time; time.sleep(5)"], timeout=1
    )
    assert returncode == 124
    assert stdout == ""
    assert "Timed out after 1s" in stderr


def test_container_smoke_builds_docker_argv():
    helper = importlib.util.spec_from_file_location(
        "container_smoke", BIN / "container_smoke.py"
    )
    assert helper is not None and helper.loader is not None
    module = importlib.util.module_from_spec(helper)
    helper.loader.exec_module(module)
    argv = module.build_argv(
        "docker", "quay.io/example/tool@sha256:" + "a" * 64,
        "tool --version", "linux/amd64", False
    )
    assert argv[:5] == ["docker", "run", "--rm", "--platform", "linux/amd64"]
    assert argv[-3:] == ["sh", "-lc", "tool --version"]


def test_bambu_manifest_records_content_checksums(tmp_path: Path):
    bam = tmp_path / "reads.bam"
    reference = tmp_path / "reference.fa"
    bam.write_bytes(b"bam")
    reference.write_text(">chr1\nAAAA\n")
    result = run_helper(
        "write_bambu_manifest.py", "--bam", str(bam), "--reference", str(reference),
        "--track-reads", "true", "--quantify", "false", "--output", str(tmp_path / "bambu.json")
    )
    assert result.returncode == 0, result.stderr
    record = json.loads((tmp_path / "bambu.json").read_text())
    assert len(record["bam"]["sha256"]) == 64
    assert len(record["reference"]["sha256"]) == 64
    assert record["track_reads"] is True
    assert record["quantify"] is False


def test_publication_checksum_verifier_detects_tampering(tmp_path: Path):
    product = tmp_path / "model.bed12"
    product.write_text("chr1\t0\t10\ttx\n")
    checksums = tmp_path / "checksums.tsv"
    import hashlib
    checksums.write_text(f"path\tsha256\n{product.name}\t{hashlib.sha256(product.read_bytes()).hexdigest()}\n")
    report = tmp_path / "verification.tsv"
    result = run_helper("verify_publication_checksums.py", "--root", str(tmp_path), "--checksums", str(checksums), "--report", str(report))
    assert result.returncode == 0, result.stderr
    product.write_text("tampered\n")
    failed = run_helper("verify_publication_checksums.py", "--root", str(tmp_path), "--checksums", str(checksums), "--report", str(report))
    assert failed.returncode != 0
    assert "ERROR" in report.read_text()


def test_publication_checksum_verifier_reports_published_path(tmp_path: Path):
    source = tmp_path / "model_reports" / "validated.bed12"
    source.parent.mkdir()
    source.write_text("chr1\t0\t10\ttx\n")
    checksums = tmp_path / "checksums.tsv"
    import hashlib
    checksums.write_text(
        "path\tsource_path\tsha256\n"
        f"per_run/run1/models/tama/normalized/run1.tama.bed12\t"
        f"model_reports/{source.name}\t{hashlib.sha256(source.read_bytes()).hexdigest()}\n"
    )
    report = tmp_path / "verification.tsv"
    result = run_helper(
        "verify_publication_checksums.py", "--root", str(tmp_path),
        "--checksums", str(checksums), "--report", str(report)
    )
    assert result.returncode == 0, result.stderr
    assert "per_run/run1/models/tama/normalized/run1.tama.bed12" in report.read_text()
    assert "model_reports/validated.bed12" in report.read_text()


def test_shard_plan_retains_empty_contig(tmp_path: Path):
    workload = tmp_path / "workload.tsv"
    workload.write_text("contig\treference_length\tmapped_reads\teligible\tresource_class\nchr1\t100\t2\ttrue\tunassigned\nchr2\t50\t0\tfalse\tunassigned\nchr3\t1000\t1000\ttrue\tunassigned\n")
    plan = tmp_path / "plan.tsv"
    result = run_helper("plan_contig_shards.py", "--workload", str(workload), "--mode", "contig", "--threshold", "10", "--very-large-threshold", "1000", "--run-id", "run1", "--plan", str(plan))
    assert result.returncode == 0, result.stderr
    text = plan.read_text()
    assert "run1__chr2" in text
    assert "EMPTY" in text
    assert "run1__chr3" in text
    assert "very_large" in text


def test_shard_coverage_allows_a_valid_zero_read_plan(tmp_path: Path):
    plan = tmp_path / "plan.tsv"
    plan.write_text("shard_id\tcontig\treference_length\tmapped_reads\teligible\tresource_class\tsource_run_id\tshard_status\nrun1__chr1\tchr1\t100\t0\tfalse\tempty\trun1\tEMPTY\n")
    report = tmp_path / "coverage.tsv"
    result = run_helper("validate_shard_coverage.py", "--plan", str(plan), "--source-count", "0", "--shard-count", "0", "--report", str(report))
    assert result.returncode == 0, result.stderr
    assert "errors\t0" in report.read_text()


def test_backend_shard_audit_reports_skipped_and_attempted_shards(tmp_path: Path):
    plan = tmp_path / "plan.tsv"
    plan.write_text("shard_id\tcontig\treference_length\tmapped_reads\teligible\tresource_class\tsource_run_id\tshard_status\nrun1__chr1\tchr1\t100\t2\ttrue\tsmall\trun1\tPLANNED\nrun1__chr2\tchr2\t50\t0\tfalse\tempty\trun1\tEMPTY\n")
    audit = tmp_path / "audit.tsv"
    audit.write_text("backend\tshard_id\tstatus\ntama\trun1__chr1\tSUCCESS\n")
    report = tmp_path / "audit.report.tsv"
    result = run_helper("audit_backend_shards.py", "--plan", str(plan), "--run-id", "run1", "--backends", "tama", "--audits", str(audit), "--report", str(report))
    assert result.returncode == 0, result.stderr
    assert "run1__chr2\tSKIPPED\tNO\tSKIPPED" in report.read_text()


def test_gtf_to_bed12_and_validation(tmp_path: Path):
    gtf = tmp_path / "models.gtf"
    gtf.write_text('chr1\tpipeline\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "tx";\nchr1\tpipeline\texon\t21\t30\t.\t+\t.\tgene_id "g"; transcript_id "tx";\n')
    bed = tmp_path / "models.bed12"
    result = run_helper("gtf_to_bed12.py", "--gtf", str(gtf), "--bed12", str(bed), "--backend", "tama", "--scope", "run")
    assert result.returncode == 0, result.stderr
    report = tmp_path / "validation.tsv"
    result = run_helper("validate_bed12.py", "--bed12", str(bed), "--report", str(report))
    assert result.returncode == 0, result.stderr
    assert "errors\t0" in report.read_text()


def test_bed12_to_gtf_preserves_strand_and_exon_blocks(tmp_path: Path):
    bed = tmp_path / "tama.bed"
    bed.write_text("chr1\t0\t30\tmodel-1\t0\t-\t0\t30\t0\t2\t10,10,\t0,20,\n")
    gtf = tmp_path / "tama.gtf"
    result = run_helper("bed12_to_gtf.py", "--bed12", str(bed), "--gtf", str(gtf), "--source", "tama")
    assert result.returncode == 0, result.stderr
    rows = [line.split("\t") for line in gtf.read_text().splitlines()]
    assert [(row[3], row[4], row[6]) for row in rows] == [("1", "10", "-"), ("21", "30", "-")]
    assert all('transcript_id "model-1"' in row[8] for row in rows)


def test_bam_junction_extraction_excludes_non_primary_reads(tmp_path: Path):
    junctions = tmp_path / "junctions.bed"
    result = run_helper("bam_to_junction_bed.py", "--bam", str(DATA / "reads.bam"), "--output", str(junctions))
    assert result.returncode == 0, result.stderr
    assert junctions.read_text().splitlines() == ["chr1\t29\t39\tchr1:29-39\t1\t+"]
    sam_junctions = tmp_path / "sam.junctions.bed"
    sam_result = run_helper("bam_to_junction_bed.py", "--sam", str(DATA / "reads.sam"), "--output", str(sam_junctions))
    assert sam_result.returncode == 0, sam_result.stderr
    assert sam_junctions.read_text().splitlines() == ["chr1\t29\t39\tchr1:29-39\t1\t+"]


def test_bam_reference_dictionary_must_match(tmp_path: Path):
    idxstats = tmp_path / "idxstats.tsv"
    idxstats.write_text("chr1\t100\t1\t0\nchr2\t50\t0\t0\n*\t0\t0\t0\n")
    dictionary = tmp_path / "reference.dict.tsv"
    dictionary.write_text("chr1\t100\nchr2\t50\n")
    report = tmp_path / "validation.tsv"
    result = run_helper("validate_bam_reference.py", "--idxstats", str(idxstats), "--reference-dictionary", str(dictionary), "--report", str(report))
    assert result.returncode == 0, result.stderr
    assert "errors\t0" in report.read_text()

    dictionary.write_text("chr1\t101\nchr2\t50\n")
    failed = run_helper("validate_bam_reference.py", "--idxstats", str(idxstats), "--reference-dictionary", str(dictionary), "--report", str(tmp_path / "bad.tsv"))
    assert failed.returncode != 0


def test_bam_to_alignment_gtf_preserves_splice_blocks_and_filters_flags(tmp_path: Path):
    gtf = tmp_path / "alignment.gtf"
    result = run_helper("bam_to_alignment_gtf.py", "--bam", str(DATA / "reads.bam"), "--output", str(gtf))
    assert result.returncode == 0, result.stderr
    rows = [line.split("\t") for line in gtf.read_text().splitlines()]
    spliced = [row for row in rows if 'transcript_id "read_spliced"' in row[8]]
    assert [(row[3], row[4]) for row in spliced] == [("20", "29"), ("40", "49")]
    assert not any('transcript_id "read_duplicate"' in row[8] for row in rows)
    sam_gtf = tmp_path / "alignment.sam.gtf"
    sam_result = run_helper("bam_to_alignment_gtf.py", "--sam", str(DATA / "reads.sam"), "--output", str(sam_gtf))
    assert sam_result.returncode == 0, sam_result.stderr
    assert 'transcript_id "read_spliced"' in sam_gtf.read_text()


def test_invalid_bed12_is_rejected(tmp_path: Path):
    bed = tmp_path / "invalid.bed12"
    bed.write_text("chr1\t0\t10\tbad\t0\t+\t0\t10\t0\t2\t10,10,\t0,0,\n")
    result = run_helper("validate_bed12.py", "--bed12", str(bed), "--report", str(tmp_path / "report.tsv"))
    assert result.returncode != 0


def test_evidence_id_changes_when_dictionary_only_policy_changes():
    payload = '{"source_dataset_id":"abc","compatibility_status":"compatible","compatibility_level":"exact","dependency_class":"model_alignment","dictionary_only":false}'
    other = payload.replace('false', 'true')
    first = run_helper("evidence_dataset_id.py", "--kind", "assessment", "--payload", payload)
    second = run_helper("evidence_dataset_id.py", "--kind", "assessment", "--payload", other)
    assert first.returncode == 0 and second.returncode == 0
    assert first.stdout.strip() != second.stdout.strip()


def test_evidence_package_builds_explicit_relation_graph(tmp_path: Path):
    models = tmp_path / "models"
    alignments = tmp_path / "alignments"
    models.mkdir()
    alignments.mkdir()
    bed12 = models / "run1.tama.run.validated.bed12"
    bed12.write_text("chr1\t0\t10\ttx\t0\t+\t0\t10\t0\t1\t10,\t0,\n")
    (alignments / "run1.alignment_manifest.json").write_text(json.dumps({"run_id": "run1", "bam_sha256": "bam"}))
    reference = tmp_path / "reference_identity.json"
    reference.write_text(json.dumps({"sha256": "reference"}))
    manifest = tmp_path / "normalized_manifest.tsv"
    manifest.write_text("id\trun_id\nrun1\trun1\n")
    records = tmp_path / "records.json"
    records.write_text(json.dumps([{
        "bed12_name": bed12.name,
        "backend": "tama",
        "scope": "run",
        "run_id": "run1",
        "cohort_id": "fixture",
        "model_manifest": "{}",
        "validation_report": "errors\\t0",
    }]))
    evidence = tmp_path / "evidence"
    packaged = run_helper(
        "package_evidence.py", "--models-dir", str(models), "--alignment-manifests-dir", str(alignments),
        "--normalized-manifest", str(manifest), "--reference-identity", str(reference),
        "--model-records", str(records), "--output-dir", str(evidence),
    )
    assert packaged.returncode == 0, packaged.stderr
    verified = run_helper("verify_evidence_publication.py", "--evidence-dir", str(evidence))
    assert verified.returncode == 0, verified.stderr
    catalog = json.loads((evidence / "catalog.json").read_text())
    assert catalog["schema"] == "long-read-evidence-v1"
    assert catalog["relations"]["rows"] == 1
    parquet_bytes = (evidence / "relations.parquet").read_bytes()
    assert parquet_bytes.startswith(b"PAR1")
    assert parquet_bytes.endswith(b"PAR1")
    assert len(list((evidence / "assessments").rglob("assessment.json"))) == 1
