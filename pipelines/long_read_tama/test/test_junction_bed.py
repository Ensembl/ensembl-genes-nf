import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "bam_to_junction_bed.py"
SPEC = importlib.util.spec_from_file_location("bam_to_junction_bed", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_junctions_ignore_clipping_and_indels():
    # 10S 30M 100N 5I 20M 2D 10M 5H
    assert MODULE.junctions_from_cigar(
        100, [(4, 10), (0, 30), (3, 100), (1, 5), (0, 20), (2, 2), (0, 10), (5, 5)]
    ) == [(130, 230)]


def test_junctions_support_multiple_introns():
    assert MODULE.junctions_from_cigar(
        0, [(0, 25), (3, 75), (0, 10), (3, 125), (0, 40)]
    ) == [(25, 100), (110, 235)]
