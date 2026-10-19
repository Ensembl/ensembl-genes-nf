from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "bin" / "report_long_read_backends.py"


def write_table(path: Path, header: str, row: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + "\n" + row + "\n")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reporting_preserves_status_and_compares_only_trace_backed_products(tmp_path: Path):
    result = tmp_path / "results"
    bed_dir = result / "native" / "cohort"
    bed_dir.mkdir(parents=True)
    bed = "chr1\t0\t100\tm1\t0\t+\t0\t100\t0\t2\t20,30\t0,70\n"
    (bed_dir / "stringtie2_combined_models.bed").write_text(bed)
    (bed_dir / "stringtie3_combined_models.bed").write_text(bed + "chr1\t200\t300\tm2\t0\t+\t200\t300\t0\t1\t100\t0\n")
    write_table(result / "reports" / "stringtie2_status.tsv", "backend\tscope\taccession_or_cohort\tstatus\tmodel_path\tmodel_count", "stringtie2\tcohort\tcohort\tSUCCESS\tnative/cohort/stringtie2_combined_models.bed\t1")
    write_table(result / "reports" / "stringtie3_status.tsv", "backend\tscope\taccession_or_cohort\tstatus\tmodel_path\tmodel_count", "stringtie3\tcohort\tcohort\tSUCCESS\tnative/cohort/stringtie3_combined_models.bed\t2")
    (result / "pipeline_info").mkdir()
    write_table(result / "pipeline_info" / "execution_trace.txt", "name\tstatus\thash\tcpus\trealtime\tpeak_rss", "STRINGTIE2_MERGE\tCOMPLETED\tabc\t2\t1m\t1 GB")
    with (result / "pipeline_info" / "execution_trace.txt").open("a") as handle:
        handle.write("STRINGTIE3_MERGE\tCOMPLETED\tdef\t2\t2m\t2 GB\n")
    out = tmp_path / "out"
    subprocess.run(["python3", str(SCRIPT), "--results-root", str(result), "--output-dir", str(out)], check=True)
    statuses = rows(out / "results" / "backend_status.tsv")
    assert {row["status"] for row in statuses if row["backend"] in {"stringtie2", "stringtie3"}} == {"SUCCESS"}
    comparisons = rows(out / "results" / "comparison_summary.tsv")
    assert len(comparisons) == 1
    assert comparisons[0]["shared_models"] == "1"
    assert comparisons[0]["denominator"] == "2"


def test_unlinked_model_is_not_reported_as_success(tmp_path: Path):
    result = tmp_path / "results"
    product = result / "published" / "tama_combined_models.bed"
    product.parent.mkdir(parents=True)
    product.write_text("chr1\t0\t10\tx\t0\t+\t0\t10\t0\t1\t10\t0\n")
    out = tmp_path / "out"
    subprocess.run(["python3", str(SCRIPT), "--results-root", str(result), "--output-dir", str(out)], check=True)
    statuses = rows(out / "results" / "backend_status.tsv")
    row = next(row for row in statuses if row["output_path"].endswith("tama_combined_models.bed"))
    assert row["status"] == "STALE_OR_UNVERIFIED"
