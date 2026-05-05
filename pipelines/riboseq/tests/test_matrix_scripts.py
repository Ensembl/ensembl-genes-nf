import gzip
import importlib.util
import pickle
import subprocess
import sys
from pathlib import Path

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


def test_hash_collision_vocab_entries_resolve_by_sequence(tmp_path):
    build_study_matrix = load_script("build_study_matrix.py")
    merge_global_matrix = load_script("merge_global_matrix.py")

    old_build_hash = build_study_matrix.hash_sequence
    old_merge_hash = merge_global_matrix.hash_sequence
    build_study_matrix.hash_sequence = lambda seq: 1
    merge_global_matrix.hash_sequence = lambda seq: 1
    try:
        tsv = tmp_path / "s1.tsv"
        tsv.write_text("AA\t2\nCC\t3\n")

        builder = build_study_matrix.StudyMatrixBuilder("STUDY", tmp_path)
        builder.build_from_tsvs(["s1"], [tsv])
        _, vocab_path, _ = builder.save()

        with open(vocab_path, "rb") as handle:
            vocab = pickle.load(handle)

        assert merge_global_matrix.resolve_local_id(vocab["seq_to_id"], "AA") == 0
        assert merge_global_matrix.resolve_local_id(vocab["seq_to_id"], "CC") == 1
    finally:
        build_study_matrix.hash_sequence = old_build_hash
        merge_global_matrix.hash_sequence = old_merge_hash


def test_merge_global_matrix_preserves_counts_across_studies(tmp_path):
    study1 = tmp_path / "study1"
    study2 = tmp_path / "study2"
    outdir = tmp_path / "global"
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

    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "merge_global_matrix.py"),
            "--study-dirs",
            str(study1),
            str(study2),
            "--output-dir",
            str(outdir),
            "--chunk-size",
            "2",
        ],
        check=True,
    )

    assert (outdir / "global_reads.fasta").read_text().splitlines() == [
        ">read_0",
        "AA",
        ">read_1",
        "CC",
        ">read_2",
        "GG",
    ]
    assert sp.load_npz(outdir / "global_matrix.npz").toarray().tolist() == [
        [2, 5],
        [3, 0],
        [0, 7],
    ]
