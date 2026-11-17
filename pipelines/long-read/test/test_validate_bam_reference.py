from __future__ import annotations

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "bin" / "validate_bam_reference.py"
SPEC = importlib.util.spec_from_file_location("validate_bam_reference", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_reordered_contigs_with_same_lengths_are_accepted(tmp_path: Path) -> None:
    reference = tmp_path / "reference.dict.tsv"
    idxstats = tmp_path / "bam.idxstats.tsv"
    report = tmp_path / "reference.validation.tsv"
    reference.write_text("chr2\t20\nchr1\t10\n")
    idxstats.write_text("chr1\t10\t0\t0\nchr2\t20\t0\t0\n*\t0\t0\t0\n")

    assert MODULE.read_dictionary(reference) != MODULE.read_idxstats(idxstats)
    assert dict(MODULE.read_dictionary(reference)) == dict(MODULE.read_idxstats(idxstats))
