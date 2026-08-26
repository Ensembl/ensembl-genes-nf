import csv
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "prepare_ribotricer_offsets.py"


def test_prepare_ribotricer_offsets(tmp_path):
    source = tmp_path / "sample.best_offset.txt"
    with source.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["length", "offset"])
        writer.writerow([30, 15])
        writer.writerow([28, 13])
    lengths = tmp_path / "read_lengths.txt"
    offsets = tmp_path / "psite_offsets.txt"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(source), "--read-lengths", str(lengths), "--offsets", str(offsets)],
        check=True,
    )
    assert lengths.read_text() == "28,30\n"
    assert offsets.read_text() == "13,15\n"
