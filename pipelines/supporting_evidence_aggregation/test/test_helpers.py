from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


def load(name: str):
    path = ROOT / "bin" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SupportingEvidenceHelpersTests(unittest.TestCase):
    def test_plan_has_one_row_per_selected_track(self) -> None:
        validator = load("validate_supporting_evidence_inputs")
        fixture_dir = ROOT.parent / "track_file_aggregation" / "test"
        plan, metadata = validator.build_plan(
            validator.read_tsv(fixture_dir / "track_manifest.tsv"),
            fixture_dir,
            {"coverage", "splice_junction"},
            "run",
            "test_species",
            False,
        )
        self.assertEqual(len(plan), 9)
        self.assertEqual(metadata["source_entity_ids"], ["run_1", "run_2", "run_3"])

    def test_splice_matrix_keeps_run_specific_support(self) -> None:
        builder = load("build_splice_matrix")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.bed"
            second = root / "second.bed"
            output = root / "matrix.tsv"
            metadata = root / "matrix.json"
            first.write_text("chr1\t10\t20\tj1\t10\t+\t3\t1\tGT-AG\n")
            second.write_text("chr1\t10\t20\tj1\t10\t+\t2\t4\tGT-AG\n")
            builder.build_matrix({"run_a": first, "run_b": second}, output, metadata)
            lines = output.read_text().splitlines()
            self.assertIn("run_a.unique_reads", lines[0])
            self.assertIn("\t3\t2\t1\t4", lines[1])

    def test_splice_matrix_rejects_conflicting_motifs(self) -> None:
        builder = load("build_splice_matrix")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.bed"
            second = root / "second.bed"
            first.write_text("chr1\t10\t20\tj1\t10\t+\t3\t1\tGT-AG\n")
            second.write_text("chr1\t10\t20\tj1\t10\t+\t2\t4\tGC-AG\n")
            with self.assertRaisesRegex(ValueError, "conflicting splice motifs"):
                builder.build_matrix({"run_a": first, "run_b": second}, root / "matrix.tsv", root / "matrix.json")


if __name__ == "__main__":
    unittest.main()
