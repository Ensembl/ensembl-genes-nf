"""Scientific contract checks for direct and model-linked support."""
from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "test" / "data"


def read_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def record(tmp_path: Path, source: str):
    support = tmp_path / f"{source}.direct_support"
    command = [
        sys.executable, str(ROOT / "bin" / "bam_to_features.py"),
        f"--{source}", str(DATA / f"reads.{source}"),
        "--read-features", str(tmp_path / f"{source}.features.tsv.gz"),
        "--five-prime", str(tmp_path / f"{source}.five.bed"),
        "--three-prime", str(tmp_path / f"{source}.three.bed"),
        "--junctions", str(tmp_path / f"{source}.junctions.bed"),
        "--summary", str(tmp_path / f"{source}.summary.tsv"),
        "--support-dir", str(support),
        "--run-id", "run1", "--assembly-id", "fixture",
        "--sample-id", "sample1", "--library-id", "lib1",
        "--technology", "pacbio", "--min-mapq", "20",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return support


def test_direct_records_match_declared_filter_and_bam_sam_agree(tmp_path: Path):
    sam = record(tmp_path, "sam")
    flow = {row["stage"]: int(row["alignment_records"]) for row in read_rows(sam / "alignment_filter_flow.tsv")}
    assert flow["records_scanned"] == 7
    assert flow["primary_mapped"] == 5
    assert flow["mapq_passing"] == 5
    assert flow["legacy_retained"] == 3
    assert flow["duplicate_records"] == flow["qc_fail_records"] == 1
    chains = read_rows(sam / "chain_observations.tsv.gz")
    assert next(row for row in chains if row["intron_chain"] == "29-39")["alignment_records"] == "1"
    ends = read_rows(sam / "end_observations.tsv.gz")
    assert next(row for row in ends if row["intron_chain"] == "29-39")["left_aligned_end"] == "19"
    assert next(row for row in ends if row["intron_chain"] == "29-39")["right_aligned_end"] == "48"
    assert next(row for row in ends if row["intron_chain"] == "" and row["alignment_orientation"] == "+")["alignment_records"] == "3"
    assert len(read_rows(tmp_path / "sam.features.tsv.gz")) == 3
    if shutil.which("samtools"):
        bam = record(tmp_path, "bam")
        for name in ("alignment_filter_flow.tsv", "chain_observations.tsv.gz", "junction_observations.tsv.gz", "end_observations.tsv.gz"):
            assert read_rows(bam / name) == read_rows(sam / name)


def bed(contig, start, end, name, strand, sizes, offsets):
    return "\t".join((contig, str(start), str(end), name, "0", strand, str(start), str(end), "0", str(len(sizes)),
                      ",".join(map(str, sizes)) + ",", ",".join(map(str, offsets)) + ",")) + "\n"


def test_model_join_preserves_exact_partial_monoexon_and_strand(tmp_path: Path):
    support = record(tmp_path, "sam")
    direct_root = tmp_path / "direct"
    direct_root.mkdir()
    support.rename(direct_root / "run1.direct_support")
    models = tmp_path / "models"
    models.mkdir()
    name = "run1.test.run.bed12"
    (models / name).write_text(
        bed("chr1", 19, 49, "test.run.model_000001.exact", "+", (10, 10), (0, 20))
        + bed("chr1", 19, 69, "test.run.model_000002.partial", "+", (10, 10, 10), (0, 20, 40))
        + bed("chr1", 4, 24, "test.run.model_000003.mono", "+", (20,), (0,))
        + bed("chr1", 54, 74, "test.run.model_000004.minus", "-", (20,), (0,))
        + bed("chr1", 54, 74, "test.run.model_000005.ambiguous", ".", (20,), (0,))
        + bed("chr1", 0, 49, "test.run.model_000006.shifted", "+", (29, 10), (0, 39))
    )
    records = tmp_path / "records.json"
    records.write_text(json.dumps([{"bed12_name": name, "backend": "test", "scope": "run", "run_id": "run1", "cohort_id": "fixture"}]))
    output = tmp_path / "joined"
    result = subprocess.run([
        sys.executable, str(ROOT / "bin" / "build_model_support.py"),
        "--models-dir", str(models), "--model-records", str(records),
        "--direct-support-dir", str(direct_root), "--reference", str(DATA / "reference.fa"),
        "--output-dir", str(output),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    index = {row["native_transcript_id"]: row["model_id"] for row in read_rows(output / "model_index.tsv.gz")}
    assert set(index) == {"exact", "partial", "mono", "minus", "ambiguous", "shifted"}
    support_rows = {row["model_id"]: row for row in read_rows(output / "model_support.tsv.gz")}
    assert support_rows[index["exact"]]["exact_chain_records"] == "1"
    assert support_rows[index["partial"]]["partial_chain_records"] == "1"
    assert support_rows[index["partial"]]["exact_chain_records"] == "0"
    assert support_rows[index["mono"]]["terminal_depth"] == "3"
    assert support_rows[index["exact"]]["read_set_id"] == "unknown"
    assert support_rows[index["exact"]]["construction_relation"] == "used_to_build"
    assert support_rows[index["ambiguous"]]["support_status"] == "ambiguous_strand"
    events = {(row["model_id"], row["event_type"]): row for row in read_rows(output / "model_event_support.tsv.gz")}
    assert (index["ambiguous"], "5prime") not in events
    assert events[index["exact"], "5prime"]["end_records_within_10bp"] == "1"
    assert events[index["minus"], "5prime"]["coordinate"] == "73"
    assert events[index["minus"], "3prime"]["coordinate"] == "54"
    assert events[index["minus"], "5prime"]["end_records_within_10bp"] == "1"
    assert events[index["exact"], "5prime"]["observed_mode_coordinate"] == "19"
    assert events[index["exact"], "5prime"]["distance_to_mode_bp"] == "0"
    assert events[index["exact"], "5prime"]["end_records_exact"] == "1"
    assert events[index["shifted"], "5prime"]["raw_chain_depth"] == "1"
    assert events[index["shifted"], "5prime"]["end_records_within_10bp"] == "0"
    assert events[index["exact"], "5prime"]["read_set_id"] == "unknown"
    junctions = [row for row in read_rows(output / "model_junction_support.tsv.gz") if row["model_id"] == index["exact"]]
    assert len(junctions) == 1 and junctions[0]["alignment_records"] == "1"
    assert len(read_rows(output / "event_support_sensitivity.tsv.gz")) > 0
    assert all("read_set_count" in row and "raw_read_set_count" in row for row in read_rows(output / "event_support.tsv.gz"))
    shifted_event = next(row for row in read_rows(output / "event_support.tsv.gz") if row["event_type"] == "5prime" and row["coordinate"] == "0")
    assert shifted_event["raw_run_count"] == "0"
    assert shifted_event["exact_run_count"] == "0"
    assert shifted_event["review_tier"] == "model_only"
    chains = {row["intron_chain"]: row for row in read_rows(output / "chain_event_support.tsv.gz")}
    assert chains["29-39"]["best_exact_chain_records"] == "1"
    assert chains["29-39"]["best_partial_chain_records"] == "0"
    assert chains["29-39;49-59"]["best_exact_chain_records"] == "0"
    assert chains["29-39;49-59"]["best_partial_chain_records"] == "1"
    assert chains["29-39;49-59"]["review_tier"] == "partial_only"
    for table_path in output.glob("*.tsv.gz"):
        with gzip.open(table_path, "rt", newline="") as handle:
            assert all(None not in row for row in csv.DictReader(handle, delimiter="\t")), table_path.name


def test_bad_mapq_is_rejected_before_recording(tmp_path: Path):
    result = subprocess.run([
        sys.executable, str(ROOT / "bin" / "bam_to_features.py"),
        "--sam", str(DATA / "reads.sam"),
        "--read-features", str(tmp_path / "features.gz"),
        "--five-prime", str(tmp_path / "five.bed"),
        "--three-prime", str(tmp_path / "three.bed"),
        "--junctions", str(tmp_path / "junctions.bed"),
        "--summary", str(tmp_path / "summary.tsv"),
        "--support-dir", str(tmp_path / "support"),
        "--run-id", "run1", "--min-mapq", "-1",
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "nonnegative" in result.stderr


def test_disabled_evidence_keeps_legacy_features_without_recording(tmp_path: Path):
    result = subprocess.run([
        sys.executable, str(ROOT / "bin" / "bam_to_features.py"),
        "--sam", str(DATA / "reads.sam"),
        "--read-features", str(tmp_path / "features.gz"),
        "--five-prime", str(tmp_path / "five.bed"),
        "--three-prime", str(tmp_path / "three.bed"),
        "--junctions", str(tmp_path / "junctions.bed"),
        "--summary", str(tmp_path / "summary.tsv"),
        "--support-dir", str(tmp_path / "direct_support"),
        "--run-id", "run1", "--skip-support",
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert len(read_rows(tmp_path / "features.gz")) == 3
    assert (tmp_path / "summary.tsv").is_file()
    assert list((tmp_path / "direct_support").iterdir()) == []


def test_corrupt_compressed_alignment_is_rejected(tmp_path: Path):
    source = tmp_path / "bad.sam.gz"
    source.write_bytes(b"not a gzip stream")
    result = subprocess.run([
        sys.executable, str(ROOT / "bin" / "bam_to_features.py"),
        "--sam", str(source),
        "--read-features", str(tmp_path / "features.gz"),
        "--five-prime", str(tmp_path / "five.bed"),
        "--three-prime", str(tmp_path / "three.bed"),
        "--junctions", str(tmp_path / "junctions.bed"),
        "--summary", str(tmp_path / "summary.tsv"),
        "--support-dir", str(tmp_path / "direct_support"),
        "--run-id", "run1",
    ], capture_output=True, text=True)
    assert result.returncode != 0
    assert "Not a gzipped file" in result.stderr


def test_support_package_detects_checksum_tampering(tmp_path: Path):
    direct = tmp_path / "direct"
    direct.mkdir()
    record(tmp_path, "sam").rename(direct / "run1.direct_support")
    models = tmp_path / "models"
    models.mkdir()
    bed_name = "run1.test.run.validated.bed12"
    (models / bed_name).write_text(bed("chr1", 19, 49, "tx", "+", (10, 10), (0, 20)))
    alignments = tmp_path / "alignments"
    alignments.mkdir()
    (alignments / "run1.alignment_manifest.json").write_text(json.dumps({"run_id": "run1", "bam_sha256": "bam"}))
    reference = tmp_path / "reference_identity.json"
    reference.write_text(json.dumps({"sha256": "reference"}))
    manifest = tmp_path / "normalized_manifest.tsv"
    manifest.write_text("id\trun_id\nrun1\trun1\n")
    records = tmp_path / "records.json"
    records.write_text(json.dumps([{"bed12_name": bed_name, "backend": "test", "scope": "run", "run_id": "run1", "cohort_id": "fixture"}]))
    joined = tmp_path / "joined"
    joined.mkdir()
    joined_names = (
        "model_index.tsv.gz", "model_support.tsv.gz", "model_junction_support.tsv.gz",
        "model_event_support.tsv.gz", "chain_event_support.tsv.gz", "event_caller_run_support.tsv.gz",
        "event_support.tsv.gz", "event_support_sensitivity.tsv.gz",
        "event_cluster_sensitivity.tsv.gz", "event_tier_summary.tsv",
        "event_best_depth_summary.tsv", "chain_length_depth_composition.tsv.gz",
        "support_schema.json",
    )
    for name in joined_names:
        (joined / name).write_text("fixture\n")
    evidence = tmp_path / "evidence"
    package = subprocess.run([
        sys.executable, str(ROOT / "bin" / "package_evidence.py"),
        "--models-dir", str(models), "--alignment-manifests-dir", str(alignments),
        "--normalized-manifest", str(manifest), "--reference-identity", str(reference),
        "--model-records", str(records), "--direct-support-dir", str(direct),
        "--joined-support-dir", str(joined), "--output-dir", str(evidence),
    ], capture_output=True, text=True)
    assert package.returncode == 0, package.stderr
    verify = [sys.executable, str(ROOT / "bin" / "verify_evidence_publication.py"),
              "--evidence-dir", str(evidence), "--require-support"]
    valid = subprocess.run(verify, capture_output=True, text=True)
    assert valid.returncode == 0, valid.stderr
    target = evidence / "support" / "summary" / "model_support.tsv.gz"
    target.write_bytes(target.read_bytes() + b"tampered")
    invalid = subprocess.run(verify, capture_output=True, text=True)
    assert invalid.returncode != 0
    assert "checksum mismatch" in invalid.stderr
