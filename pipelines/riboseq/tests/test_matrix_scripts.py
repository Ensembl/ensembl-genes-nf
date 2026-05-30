import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import polars as pl
import scipy.sparse as sp


RIBOSEQ_DIR = Path(__file__).resolve().parents[1]
BIN_DIR = RIBOSEQ_DIR / "bin"


def load_script(name):
    path = BIN_DIR / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_collapsed_to_tsv_spills_sorts_and_coalesces_duplicates(tmp_path):
    fasta = tmp_path / "sample.collapsed.fa"
    fasta.write_text(
        ">read_x2\n"
        "TT\n"
        ">read_x3\n"
        "AA\n"
        ">read_x4\n"
        "AA\n"
        ">missing_count\n"
        "CC\n"
        ">multi_x5\n"
        "A\n"
        "C\n"
    )
    output = tmp_path / "sample.tsv"

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "collapsed_to_tsv.py"),
            str(fasta),
            "--output",
            str(output),
            "--sample-id",
            "sample",
            "--chunk-size",
            "2",
            "--merge-fan-in",
            "2",
        ],
        check=True,
    )

    assert output.read_text().splitlines() == [
        "AA\t7",
        "AC\t5",
        "CC\t1",
        "TT\t2",
    ]


def test_collapsed_to_tsv_supports_longer_prefix_partitions(tmp_path):
    fasta = tmp_path / "sample.collapsed.fa"
    fasta.write_text(
        ">read_x2\n"
        "AAAA\n"
        ">read_x3\n"
        "AAAC\n"
        ">read_x5\n"
        "ANNN\n"
    )

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "collapsed_to_tsv.py"),
            str(fasta),
            "--partition",
            "--partition-prefix-length",
            "3",
            "--output-dir",
            str(tmp_path),
            "--sample-id",
            "sample",
            "--chunk-size",
            "65",
        ],
        check=True,
    )

    assert (tmp_path / "sample.AAA.tsv").read_text().splitlines() == [
        "AAAA\t2",
        "AAAC\t3",
    ]
    assert (tmp_path / "sample.NNN.tsv").read_text().splitlines() == ["ANNN\t5"]

    stats = json.loads((tmp_path / "sample.partition_stats.json").read_text())
    assert stats["sample_id"] == "sample"
    assert stats["partition_prefix_length"] == 3
    assert stats["catch_all_partition"] == "NNN"
    assert stats["total_records"] == 3
    assert stats["total_counts"] == 10
    assert stats["total_unique_sequences"] == 3
    partition_stats = {item["partition"]: item for item in stats["partitions"]}
    assert partition_stats["AAA"] == {
        "partition": "AAA",
        "records": 2,
        "counts": 5,
        "unique_sequences": 2,
    }
    assert partition_stats["NNN"] == {
        "partition": "NNN",
        "records": 1,
        "counts": 5,
        "unique_sequences": 1,
    }


def test_qc_partitioned_tsv_fails_empty_sample(tmp_path):
    stats = tmp_path / "sample.partition_stats.json"
    report = tmp_path / "sample.partition_qc.json"
    stats.write_text(
        json.dumps(
            {
                "sample_id": "sample",
                "partition_prefix_length": 4,
                "catch_all_partition": "NNNN",
                "total_records": 0,
                "total_counts": 0,
                "total_unique_sequences": 0,
                "partitions": [],
            }
        )
    )

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "qc_partitioned_tsv.py"),
            str(stats),
            "--output",
            str(report),
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    qc = json.loads(report.read_text())
    assert result.returncode == 1
    assert qc["passed"] is False
    assert "total_records=0 is below 1" in qc["failures"]
    assert "total_counts=0 is below 1" in qc["failures"]


def test_qc_partitioned_tsv_flags_high_catch_all_fraction(tmp_path):
    stats = tmp_path / "sample.partition_stats.json"
    report = tmp_path / "sample.partition_qc.json"
    stats.write_text(
        json.dumps(
            {
                "sample_id": "sample",
                "partition_prefix_length": 4,
                "catch_all_partition": "NNNN",
                "total_records": 10,
                "total_counts": 100,
                "total_unique_sequences": 10,
                "partitions": [
                    {"partition": "AAAA", "records": 9, "counts": 90, "unique_sequences": 9},
                    {"partition": "NNNN", "records": 1, "counts": 10, "unique_sequences": 1},
                ],
            }
        )
    )

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "qc_partitioned_tsv.py"),
            str(stats),
            "--output",
            str(report),
            "--max-catch-all-unique-fraction",
            "0.05",
            "--max-catch-all-count-fraction",
            "0.05",
        ],
        check=False,
    )

    qc = json.loads(report.read_text())
    assert result.returncode == 1
    assert qc["passed"] is False
    assert qc["observed"]["catch_all_unique_fraction"] == 0.1
    assert qc["observed"]["catch_all_count_fraction"] == 0.1


def test_qc_partitioned_tsv_record_action_keeps_report_but_exits_zero(tmp_path):
    stats = tmp_path / "sample.partition_stats.json"
    report = tmp_path / "sample.partition_qc.json"
    stats.write_text(
        json.dumps(
            {
                "sample_id": "sample",
                "partition_prefix_length": 4,
                "catch_all_partition": "NNNN",
                "total_records": 10,
                "total_counts": 100,
                "total_unique_sequences": 10,
                "partitions": [
                    {"partition": "AAAA", "records": 9, "counts": 90, "unique_sequences": 9},
                    {"partition": "NNNN", "records": 1, "counts": 10, "unique_sequences": 1},
                ],
            }
        )
    )

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "qc_partitioned_tsv.py"),
            str(stats),
            "--output",
            str(report),
            "--max-catch-all-unique-fraction",
            "0.05",
            "--max-catch-all-count-fraction",
            "0.05",
            "--action",
            "record",
        ],
        check=False,
        text=True,
        capture_output=True,
    )

    qc = json.loads(report.read_text())
    assert result.returncode == 0
    assert qc["passed"] is False
    assert qc["action"] == "record"
    assert "NNNN count fraction 0.100000 exceeds 0.050000" in qc["failures"]
    assert "ERROR: sample: NNNN count fraction 0.100000 exceeds 0.050000" in result.stderr


def test_build_study_matrix_preserves_sample_columns_and_coalesces_rows(tmp_path):
    s1 = tmp_path / "s1.tsv"
    s2 = tmp_path / "s2.tsv"
    s1.write_text("AA\t2\nCC\t3\nCC\t4\nTT\t5\n")
    s2.write_text("AA\t11\nGG\t13\nTT\t17\n")

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_study_matrix.py"),
            str(s1),
            str(s2),
            "--output-dir",
            str(tmp_path),
            "--study-id",
            "STUDY",
            "--sample-ids",
            "s1,s2",
        ],
        check=True,
    )

    with gzip.open(tmp_path / "STUDY_sequences.txt.gz", "rt") as handle:
        sequences = [line.strip() for line in handle if line.strip()]
    matrix = sp.load_npz(tmp_path / "STUDY_matrix.npz").toarray().tolist()

    assert sequences == ["AA", "CC", "GG", "TT"]
    assert matrix == [
        [2, 11],
        [7, 0],
        [0, 13],
        [5, 17],
    ]


def test_build_study_matrix_writes_empty_partition(tmp_path):
    s1 = tmp_path / "s1.NNNN.tsv"
    s2 = tmp_path / "s2.NNNN.tsv"
    s1.write_text("")
    s2.write_text("")

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_study_matrix.py"),
            str(s1),
            str(s2),
            "--output-dir",
            str(tmp_path),
            "--study-id",
            "STUDY",
            "--sample-ids",
            "s1,s2",
            "--partition",
            "NNNN",
        ],
        check=True,
    )

    with gzip.open(tmp_path / "STUDY.NNNN_sequences.txt.gz", "rt") as handle:
        assert handle.read() == ""

    matrix = sp.load_npz(tmp_path / "STUDY.NNNN_matrix.npz")
    metadata = json.loads((tmp_path / "STUDY.NNNN_metadata.json").read_text())

    assert matrix.shape == (0, 2)
    assert matrix.nnz == 0
    assert metadata["n_reads"] == 0
    assert metadata["nnz"] == 0


def make_two_study_fixture(tmp_path):
    study1 = tmp_path / "study1"
    study2 = tmp_path / "study2"
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

    for vocab in list(study1.glob("*_vocab.pkl")) + list(study2.glob("*_vocab.pkl")):
        vocab.unlink()

    return study1, study2


def test_merge_global_matrix_preserves_counts_across_studies(tmp_path):
    study1, study2 = make_two_study_fixture(tmp_path)
    outdir = tmp_path / "global"

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir),
            "--matrix-format",
            "dense",
            "--chunk-size",
            "2",
            "--metadata-shard-rows",
            "2",
        ],
        check=True,
    )

    metadata_dir = outdir / "global_metadata.parquet"
    assert metadata_dir.is_dir()
    assert sorted(path.name for path in metadata_dir.glob("part-*.parquet")) == [
        "part-00000.parquet",
        "part-00001.parquet",
    ]
    assert (outdir / "global_reads.fasta").read_text().splitlines() == [
        ">read_0",
        "AA",
        ">read_1",
        "CC",
        ">read_2",
        "GG",
    ]
    if (outdir / "global_matrix.npz").exists():
        dense_matrix = sp.load_npz(outdir / "global_matrix.npz").toarray().tolist()
    else:
        import zarr

        dense_matrix = zarr.open_array(str(outdir / "global_matrix.zarr" / "counts"), mode="r")[:].tolist()
    assert dense_matrix == [[2, 5], [3, 0], [0, 7]]

    outdir_no_fasta = tmp_path / "global_no_fasta"
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir_no_fasta),
            "--matrix-format",
            "dense",
            "--chunk-size",
            "2",
            "--metadata-shard-rows",
            "2",
            "--no-write-fasta",
        ],
        check=True,
    )

    assert (outdir_no_fasta / "global_reads.fasta").read_text() == ""
    config = json.loads((outdir_no_fasta / "global_config.json").read_text())
    assert config["write_fasta"] is False


def test_sparse_parquet_global_store_preserves_counts_across_studies(tmp_path):
    study1, study2 = make_two_study_fixture(tmp_path)
    outdir = tmp_path / "global_sparse"

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir),
            "--sparse-shard-rows",
            "2",
            "--sparse-read-bucket-size",
            "2",
            "--metadata-shard-rows",
            "2",
        ],
        check=True,
    )

    facts = (
        pl.read_parquet(str(outdir / "global_matrix.parquet" / "read_bucket=*" / "*.parquet"))
        .select("read_id", "sample_id", "study_id_int", "count")
        .sort("read_id", "sample_id")
        .to_dicts()
    )
    assert facts == [
        {"read_id": 0, "sample_id": 0, "study_id_int": 0, "count": 2},
        {"read_id": 0, "sample_id": 1, "study_id_int": 1, "count": 5},
        {"read_id": 1, "sample_id": 0, "study_id_int": 0, "count": 3},
        {"read_id": 2, "sample_id": 1, "study_id_int": 1, "count": 7},
    ]

    samples = pl.read_parquet(outdir / "global_samples.parquet").sort("sample_id").to_dicts()
    assert samples == [
        {
            "sample_id": 0,
            "sample_name": "s1",
            "study_id_int": 0,
            "study_id": "STUDY1",
            "study_sample_index": 0,
        },
        {
            "sample_id": 1,
            "sample_name": "s2",
            "study_id_int": 1,
            "study_id": "STUDY2",
            "study_sample_index": 0,
        },
    ]
    studies = pl.read_parquet(outdir / "global_studies.parquet").sort("study_id_int").to_dicts()
    assert studies == [
        {"study_id_int": 0, "study_id": "STUDY1", "n_reads": 2, "sample_offset": 0},
        {"study_id_int": 1, "study_id": "STUDY2", "n_reads": 2, "sample_offset": 1},
    ]

    manifest = json.loads((outdir / "global_matrix_manifest.json").read_text())
    assert manifest["matrix_format"] == "sparse-parquet"
    assert manifest["status"] == "committed"
    assert manifest["generation"] == 1
    assert manifest["parents"] == []
    assert manifest["tombstones_path"] == "global_tombstones.parquet"
    assert manifest["global_counts"] == "global_counts"
    assert manifest["global_reads"] == "global_reads.parquet"
    assert manifest["append_protocol"]["mode"] == "generation"
    assert manifest["schema"]["sample_id"] == "uint32"
    assert manifest["schema"]["study_id_int"] == "uint32"
    assert manifest["lookup_tables"] == {
        "samples": "global_samples.parquet",
        "studies": "global_studies.parquet",
    }
    assert manifest["n_reads"] == 3
    assert manifest["n_samples"] == 2
    assert manifest["nnz"] == 4
    assert len(manifest["generations"]) == 1
    assert manifest["generations"][0]["path"] == "global_counts/generation=000001"
    assert manifest["partitioning"] == ["generation", "read_bucket"]
    assert manifest["read_bucket_size"] == 2
    assert sorted({part["read_bucket"] for part in manifest["parts"]}) == [0, 1]
    assert sum(part["rows"] for part in manifest["parts"]) == 4

    config = json.loads((outdir / "global_config.json").read_text())
    assert config["matrix_format"] == "sparse-parquet"
    assert config["matrix_manifest"] == "global_matrix_manifest.json"
    assert config["retained_counts_materialized"] is False
    assert (outdir / "global_reads.parquet").exists()
    assert pl.read_parquet(outdir / "global_reads.parquet").sort("read_id").to_dicts() == [
        {"read_id": 0, "read_bucket": 0, "sequence": "AA", "length": 2},
        {"read_id": 1, "read_bucket": 0, "sequence": "CC", "length": 2},
        {"read_id": 2, "read_bucket": 1, "sequence": "GG", "length": 2},
    ]
    assert (outdir / "global_counts" / "generation=000001").exists()
    assert (outdir / "global_matrix.parquet").is_symlink()
    assert not (outdir / "global_retained_counts.parquet").exists()


def test_sparse_parquet_global_store_does_not_materialize_study_csr(tmp_path, monkeypatch):
    script = load_script("merge_global_matrix.py")
    study1, study2 = make_two_study_fixture(tmp_path)
    outdir = tmp_path / "global_sparse"

    def fail_load_npz(_path):
        raise AssertionError("sparse merge should stream CSR npz rows")

    monkeypatch.setattr(script.sp, "load_npz", fail_load_npz)

    study_inputs = []
    study_matrices = {}
    for study_dir in [study1, study2]:
        metadata_path = next(study_dir.glob("*_metadata.json"))
        study_id = metadata_path.stem.removesuffix("_metadata")
        n_reads, samples = script.load_study_metadata(metadata_path, study_id)
        matrix_path = study_dir / f"{study_id}_matrix.npz"
        study_inputs.append(script.StudyInput(
            study_id=study_id,
            sequences_path=study_dir / f"{study_id}_sequences.txt.gz",
            matrix_path=matrix_path,
            metadata_path=metadata_path,
            n_reads=n_reads,
            samples=samples,
        ))
        study_matrices[study_id] = matrix_path

    merger = script.GlobalMatrixMerger(
        outdir,
        matrix_format="sparse-parquet",
        sparse_read_bucket_size=2,
    )
    with script.make_tempdir("global_merge_", outdir) as temp_dir:
        try:
            merger.build_global_vocabulary(study_inputs, Path(temp_dir))
            merger.build_global_matrix(study_matrices)
        finally:
            merger.close_remaps()
    merger.save_outputs()
    assert merger.new_read_rows == []
    assert merger.existing_reads == {}

    facts = (
        pl.read_parquet(str(outdir / "global_counts" / "generation=000001" / "read_bucket=*" / "*.parquet"))
        .select("read_id", "sample_id", "count")
        .sort("read_id", "sample_id")
        .to_dicts()
    )
    assert facts == [
        {"read_id": 0, "sample_id": 0, "count": 2},
        {"read_id": 0, "sample_id": 1, "count": 5},
        {"read_id": 1, "sample_id": 0, "count": 3},
        {"read_id": 2, "sample_id": 1, "count": 7},
    ]


def test_csr_npz_streaming_uses_public_numpy_header_reader(tmp_path, monkeypatch):
    script = load_script("merge_global_matrix.py")
    study1, _ = make_two_study_fixture(tmp_path)
    matrix_path = study1 / "STUDY1_matrix.npz"

    if hasattr(script.np.lib.format, "_read_array_header"):
        monkeypatch.delattr(script.np.lib.format, "_read_array_header")

    rows = [
        (local_row, sample_indices.tolist(), counts.tolist())
        for local_row, sample_indices, counts in script.iter_csr_npz_rows(matrix_path)
    ]

    assert rows == [
        (0, [0], [2]),
        (1, [0], [3]),
    ]


def build_partition_fixture(tmp_path, partition, rows1, rows2):
    study1 = tmp_path / f"study1_{partition}"
    study2 = tmp_path / f"study2_{partition}"
    study1.mkdir()
    study2.mkdir()

    s1 = tmp_path / f"s1.{partition}.tsv"
    s2 = tmp_path / f"s2.{partition}.tsv"
    s1.write_text(rows1)
    s2.write_text(rows2)

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
            "--partition",
            partition,
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
            "--partition",
            partition,
        ],
        check=True,
    )
    return study1, study2


def test_partitioned_sparse_merge_assigns_deterministic_noncolliding_read_ids(tmp_path):
    aa_studies = build_partition_fixture(tmp_path, "AA", "AAAA\t2\nAAAC\t3\n", "AAAC\t5\n")
    cc_studies = build_partition_fixture(tmp_path, "CC", "CCCC\t7\n", "CCCG\t11\n")

    aa_out = tmp_path / "global_AA"
    cc_out = tmp_path / "global_CC"
    for outdir, studies, partition, ordinal in [
        (aa_out, aa_studies, "AA", 0),
        (cc_out, cc_studies, "CC", 1),
    ]:
        subprocess.run(
            [
                sys.executable,
                str(BIN_DIR / "merge_global_matrix.py"),
                "--study-dirs",
                str(studies[0]),
                str(studies[1]),
                "--output-dir",
                str(outdir),
                "--partition",
                partition,
                "--partition-ordinal",
                str(ordinal),
                "--partition-stride",
                "1000",
                "--sparse-read-bucket-size",
                "100",
                "--metadata-shard-rows",
                "2",
            ],
            check=True,
        )

    aa_reads = pl.read_parquet(aa_out / "global.AA_reads.parquet").sort("read_id").to_dicts()
    cc_reads = pl.read_parquet(cc_out / "global.CC_reads.parquet").sort("read_id").to_dicts()
    assert aa_reads == [
        {"read_id": 0, "read_bucket": 0, "sequence": "AAAA", "length": 4},
        {"read_id": 1, "read_bucket": 0, "sequence": "AAAC", "length": 4},
    ]
    assert cc_reads == [
        {"read_id": 1000, "read_bucket": 10, "sequence": "CCCC", "length": 4},
        {"read_id": 1001, "read_bucket": 10, "sequence": "CCCG", "length": 4},
    ]

    aa_manifest = json.loads((aa_out / "global.AA_matrix_manifest.json").read_text())
    cc_manifest = json.loads((cc_out / "global.CC_matrix_manifest.json").read_text())
    assert aa_manifest["read_id_offset"] == 0
    assert aa_manifest["partition_ordinal"] == 0
    assert aa_manifest["partition_stride"] == 1000
    assert aa_manifest["n_reads"] == 2
    assert cc_manifest["read_id_offset"] == 1000
    assert cc_manifest["partition_ordinal"] == 1
    assert cc_manifest["partition_stride"] == 1000
    assert cc_manifest["n_reads"] == 2


def test_sparse_parquet_append_writes_new_generation_with_stable_ids(tmp_path):
    study1, study2 = make_two_study_fixture(tmp_path)
    outdir = tmp_path / "global_sparse"
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
        ],
        check=True,
    )

    study3 = tmp_path / "study3"
    study3.mkdir()
    s3 = tmp_path / "s3.tsv"
    s3.write_text("AA\t11\nTT\t13\n")
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "build_study_matrix.py"),
            str(s3),
            "--output-dir",
            str(study3),
            "--study-id",
            "STUDY3",
            "--sample-ids",
            "s3",
        ],
        check=True,
    )

    append_out = tmp_path / "global_sparse_append"
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study3),
            "--output-dir",
            str(append_out),
            "--append-to",
            str(outdir),
            "--sparse-read-bucket-size",
            "2",
        ],
        check=True,
    )

    manifest = json.loads((append_out / "global_matrix_manifest.json").read_text())
    assert manifest["generation"] == 2
    assert manifest["parents"] == [1]
    assert [item["generation"] for item in manifest["generations"]] == [1, 2]
    assert (append_out / "global_counts" / "generation=000002").exists()

    reads = pl.read_parquet(append_out / "global_reads.parquet").sort("read_id").to_dicts()
    assert reads == [
        {"read_id": 0, "read_bucket": 0, "sequence": "AA", "length": 2},
        {"read_id": 1, "read_bucket": 0, "sequence": "CC", "length": 2},
        {"read_id": 2, "read_bucket": 1, "sequence": "GG", "length": 2},
        {"read_id": 3, "read_bucket": 1, "sequence": "TT", "length": 2},
    ]

    samples = pl.read_parquet(append_out / "global_samples.parquet").sort("sample_id").to_dicts()
    assert [row["sample_name"] for row in samples] == ["s1", "s2", "s3"]

    facts = (
        pl.concat([
            pl.read_parquet(str(outdir / item["path"] / "read_bucket=*" / "*.parquet"))
            if item["generation"] == 1
            else pl.read_parquet(str(append_out / item["path"] / "read_bucket=*" / "*.parquet"))
            for item in manifest["generations"]
        ])
        .select("read_id", "sample_id", "count")
        .sort("read_id", "sample_id")
        .to_dicts()
    )
    assert facts == [
        {"read_id": 0, "sample_id": 0, "count": 2},
        {"read_id": 0, "sample_id": 1, "count": 5},
        {"read_id": 0, "sample_id": 2, "count": 11},
        {"read_id": 1, "sample_id": 0, "count": 3},
        {"read_id": 2, "sample_id": 1, "count": 7},
        {"read_id": 3, "sample_id": 2, "count": 13},
    ]

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "query_sparse_global_matrix.py"),
            "query",
            str(append_out),
            "--read-ids",
            "0,3",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    assert result.stdout.splitlines() == [
        "read_id\tsample_index\tsample_id\tstudy_id\tcount",
        "0\t0\ts1\tSTUDY1\t2",
        "0\t1\ts2\tSTUDY2\t5",
        "0\t2\ts3\tSTUDY3\t11",
        "3\t2\ts3\tSTUDY3\t13",
    ]


def test_dense_global_matrix_guard_fails_before_writing_unsafe_chunk(tmp_path):
    study1, study2 = make_two_study_fixture(tmp_path)
    outdir = tmp_path / "global_dense_guard"

    result = subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir),
            "--matrix-format",
            "dense",
            "--chunk-size",
            "2",
            "--dense-chunk-byte-limit",
            "1",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "Dense global matrix output is unsafe" in result.stderr
    assert not (outdir / "global_matrix.npz").exists()
    assert not (outdir / "global_matrix.zarr").exists()
