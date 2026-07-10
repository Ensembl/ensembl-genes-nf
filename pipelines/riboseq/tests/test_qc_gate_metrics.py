import csv
import json
import subprocess
import sys
from pathlib import Path


RIBOSEQ_DIR = Path(__file__).resolve().parents[1]
BIN_DIR = RIBOSEQ_DIR / "bin"
FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "qc"
RULES = RIBOSEQ_DIR / "resources" / "qc_rules.default.yaml"


def run_collect(tmp_path, sample, json_name, csv_name):
    prefix = tmp_path / sample
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "collect_qc_metrics.py"),
            "--run-id",
            "run1",
            "--sample-id",
            sample,
            "--star-log",
            str(FIXTURE_DIR / "Log.final.out"),
            "--getrpf-report",
            str(FIXTURE_DIR / "getrpf_report.txt"),
            "--getrpf-checks",
            str(FIXTURE_DIR / "getrpf_checks.tsv"),
            "--ribometric-json",
            str(FIXTURE_DIR / json_name),
            "--ribometric-csv",
            str(FIXTURE_DIR / csv_name),
            "--offsets",
            str(FIXTURE_DIR / "offsets.tsv"),
            "--out-prefix",
            str(prefix),
        ],
        check=True,
    )
    return prefix.with_suffix(".qc_metrics.tsv")


def read_metrics(path):
    with open(path, newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {
            (row["metric"], row["scope"], row["length"]): float(row["value_num"])
            for row in rows
        }


def run_gate(tmp_path, sample, metrics_tsv):
    prefix = tmp_path / f"{sample}.gate"
    subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "qc_gate.py"),
            "--run-id",
            "run1",
            "--sample-id",
            sample,
            "--metrics-tsv",
            str(metrics_tsv),
            "--best-offset",
            str(FIXTURE_DIR / "offsets.tsv"),
            "--rules-yaml",
            str(RULES),
            "--out-prefix",
            str(prefix),
        ],
        check=True,
    )
    return prefix


def test_collect_qc_metrics_emits_ribometric_metrics_from_json(tmp_path):
    metrics_tsv = run_collect(tmp_path, "strong", "strong_RiboMetric.json", "strong_RiboMetric.csv")
    metrics = read_metrics(metrics_tsv)

    assert metrics[("ribometric.information_content", "sample", "")] == 0.50
    assert metrics[("ribometric.reads_at_length", "length", "29")] == 10000
    assert metrics[("ribometric.frame0_frac", "length", "29")] == 0.72
    assert metrics[("ribometric.periodicity_dominance", "length", "29")] == 0.72
    assert metrics[("ribometric.n_recommended_read_lengths", "sample", "")] == 1
    assert round(metrics[("ribometric.recommended_read_proportion", "sample", "")], 4) == 0.5556
    assert metrics[("ribometric.recommended_length", "length", "29")] == 1
    assert metrics[("ribometric.recommended_length", "length", "30")] == 0


def test_collect_qc_metrics_accepts_csv_read_len_alias_and_frame_counts(tmp_path):
    metrics_tsv = run_collect(
        tmp_path,
        "csv_alias",
        "csv_read_len_RiboMetric.json",
        "csv_read_len_RiboMetric.csv",
    )
    metrics = read_metrics(metrics_tsv)

    assert metrics[("ribometric.information_content", "sample", "")] == 0.45
    assert metrics[("ribometric.reads_at_length", "length", "28")] == 6000
    assert metrics[("ribometric.frame0_frac", "length", "28")] == 0.75


def test_qc_gate_selects_only_offsets_for_passing_lengths(tmp_path):
    metrics_tsv = run_collect(tmp_path, "strong", "strong_RiboMetric.json", "strong_RiboMetric.csv")
    gate_prefix = run_gate(tmp_path, "strong", metrics_tsv)

    assert (tmp_path / "strong.gate.offsets.selected.tsv").read_text().splitlines() == [
        "read_len\toffset",
        "29\t12",
    ]
    qc = json.loads((tmp_path / "strong.gate.qc.json").read_text())
    assert qc["selected_for_translon"] is True
    assert qc["pass_lengths"] == [29]


def test_qc_gate_emits_no_selected_offsets_when_periodicity_fails(tmp_path):
    metrics_tsv = run_collect(
        tmp_path,
        "SRR11005904",
        "weak_SRR11005904_RiboMetric.json",
        "weak_SRR11005904_RiboMetric.csv",
    )
    gate_prefix = run_gate(tmp_path, "SRR11005904", metrics_tsv)

    assert not Path(f"{gate_prefix}.offsets.selected.tsv").exists()
    assert Path(f"{gate_prefix}.offsets.pass.tsv").read_text().splitlines() == ["read_len\toffset"]
    qc = json.loads(Path(f"{gate_prefix}.qc.json").read_text())
    assert qc["selected_for_translon"] is False
    assert qc["pass_lengths"] == []


def test_qc_gate_empty_ribometric_report_has_no_pass_lengths(tmp_path):
    metrics_tsv = run_collect(tmp_path, "empty", "empty_RiboMetric.json", "empty_RiboMetric.csv")
    gate_prefix = run_gate(tmp_path, "empty", metrics_tsv)

    assert not Path(f"{gate_prefix}.offsets.selected.tsv").exists()
    assert Path(f"{gate_prefix}.pass_lengths.tsv").read_text().splitlines() == ["length"]
