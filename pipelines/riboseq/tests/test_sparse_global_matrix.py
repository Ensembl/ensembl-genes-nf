import json
import importlib.util
import subprocess
import sys
from pathlib import Path

import polars as pl
import pysam


RIBOSEQ_DIR = Path(__file__).resolve().parents[1]
BIN_DIR = RIBOSEQ_DIR / "bin"


def load_script(name):
    path = BIN_DIR / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_sparse_store(tmp_path):
    study1 = tmp_path / "study1"
    study2 = tmp_path / "study2"
    outdir = tmp_path / "global_sparse"
    study1.mkdir()
    study2.mkdir()

    s1 = tmp_path / "s1.tsv"
    s2 = tmp_path / "s2.tsv"
    s1.write_text("AA\t2\nCC\t3\n")
    s2.write_text("AA\t5\nGG\t7\n")

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_study_matrix.py"),
            str(s1),
            "--output-dir",
            str(study1),
            "--study-id",
            "STUDY1",
            "--sample-ids",
            "s1",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_study_matrix.py"),
            str(s2),
            "--output-dir",
            str(study2),
            "--study-id",
            "STUDY2",
            "--sample-ids",
            "s2",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir),
            "--sparse-read-bucket-size",
            "2",
            "--sparse-shard-rows",
            "10",
        ],
        check=True,
    )

    return outdir


def write_test_bam(path: Path):
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for query_name, start, flag in [
            ("read_0", 10, 0),
            ("read_1", 10, 0),
            ("read_2", 12, 16),
            ("unparseable", 20, 0),
        ]:
            record = pysam.AlignedSegment()
            record.query_name = query_name
            record.query_sequence = "A" * 28
            record.flag = flag
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 255
            record.cigartuples = [(0, 28)]
            record.query_qualities = pysam.qualitystring_to_array("I" * 28)
            bam.write(record)


def test_sparse_global_matrix_query_filters_by_read_and_sample(tmp_path):
    store = make_sparse_store(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "query",
            str(store),
            "--read-ids",
            "0,2",
            "--sample-ids",
            "s2",
        ],
        check=True,
        text=True,
        capture_output=True,
    )

    assert result.stdout.splitlines() == [
        "read_id\tsample_index\tsample_id\tstudy_id\tcount",
        "0\t1\ts2\tSTUDY2\t5",
        "2\t1\ts2\tSTUDY2\t7",
    ]


def test_sparse_global_matrix_tombstone_is_logical_and_query_respects_it(tmp_path):
    store = make_sparse_store(tmp_path)

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "tombstone",
            str(store),
            "--sample-id",
            "s2",
            "--reason",
            "test removal",
        ],
        check=True,
    )

    manifest = json.loads((store / "global_matrix_manifest.json").read_text())
    assert manifest["tombstones_path"] == "global_tombstones.parquet"
    tombstones = pl.read_parquet(store / "global_tombstones.parquet").to_dicts()
    assert tombstones[0]["sample_id"] == "s2"
    assert tombstones[0]["active"] is True
    assert (store / "global_matrix.parquet").exists()

    filtered = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "query",
            str(store),
            "--read-ids",
            "0",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert filtered.stdout.splitlines() == [
        "read_id\tsample_index\tsample_id\tstudy_id\tcount",
        "0\t0\ts1\tSTUDY1\t2",
    ]

    unfiltered = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "query",
            str(store),
            "--read-ids",
            "0",
            "--include-tombstoned",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert unfiltered.stdout.splitlines() == [
        "read_id\tsample_index\tsample_id\tstudy_id\tcount",
        "0\t0\ts1\tSTUDY1\t2",
        "0\t1\ts2\tSTUDY2\t5",
    ]


def test_sparse_global_matrix_retained_view_materializes_without_tombstoned_rows(tmp_path):
    store = make_sparse_store(tmp_path)
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "tombstone",
            str(store),
            "--sample-id",
            "s2",
            "--reason",
            "test removal",
        ],
        check=True,
    )
    out = tmp_path / "retained.parquet"
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "retained-view",
            str(store),
            "--output",
            str(out),
        ],
        check=True,
    )

    rows = pl.read_parquet(out).sort("read_id").to_dicts()
    assert rows == [
        {"read_id": 0, "sample_id": 0, "study_id": "STUDY1", "count": 2},
        {"read_id": 1, "sample_id": 0, "study_id": "STUDY1", "count": 3},
    ]


def test_profile_from_global_bam_joins_counts_by_bucket(tmp_path, monkeypatch):
    store = make_sparse_store(tmp_path)
    bam = tmp_path / "global.bam"
    write_test_bam(bam)

    script = load_script("profile_from_global_bam.py")
    load_calls = []
    original_load_count_bucket = script.load_count_bucket

    def counting_load_count_bucket(generation_roots, read_bucket):
        load_calls.append(read_bucket)
        return original_load_count_bucket(generation_roots, read_bucket)

    monkeypatch.setattr(script, "load_count_bucket", counting_load_count_bucket)

    out = tmp_path / "profile.parquet"
    stats = script.build_profile(
        bam=bam,
        manifest=store / "global_matrix_manifest.json",
        output=out,
        keep_bucket_shards=False,
        no_final_coalesce=False,
    )

    assert load_calls == [0, 1]
    assert stats["bam_events"] == 3
    assert stats["skipped_unparsed"] == 1
    assert stats["touched_buckets"] == 2
    assert stats["loaded_count_buckets"] == 2

    rows = (
        pl.read_parquet(out)
        .sort("sample_id", "chrom", "position", "strand")
        .to_dicts()
    )
    assert rows == [
        {
            "sample_id": 0,
            "study_id_int": 0,
            "chrom": "chr1",
            "position": 10,
            "strand": "+",
            "count": 5,
        },
        {
            "sample_id": 1,
            "study_id_int": 1,
            "chrom": "chr1",
            "position": 10,
            "strand": "+",
            "count": 5,
        },
        {
            "sample_id": 1,
            "study_id_int": 1,
            "chrom": "chr1",
            "position": 12,
            "strand": "-",
            "count": 7,
        },
    ]
