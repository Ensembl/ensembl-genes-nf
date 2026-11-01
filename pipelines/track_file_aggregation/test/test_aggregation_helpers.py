#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).parents[1]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


AGGREGATOR = load("aggregate_splice_junctions", ROOT / "bin/aggregate_splice_junctions.py")
VALIDATOR = load("validate_aggregation_inputs", ROOT / "bin/validate_aggregation_inputs.py")


class AggregationHelpersTests(unittest.TestCase):
    def test_sums_support_and_caps_display_score(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.bed"
            second = root / "second.bed"
            output = root / "aggregate.bed"
            first.write_text("chr1\t10\t20\tj1\t1000\t+\t900\t2\tGT-AG\n")
            second.write_text("chr1\t10\t20\tj1\t1000\t+\t200\t3\tGT-AG\n")
            AGGREGATOR.aggregate([first, second], output)
            fields = output.read_text().strip().split("\t")
            self.assertEqual(fields[6:9], ["1100", "5", "GT-AG"])
            self.assertEqual(fields[4], "1000")

    def test_conflicting_motifs_are_detectable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.bed"
            second = root / "second.bed"
            output = root / "aggregate.bed"
            first.write_text("chr1\t10\t20\tj1\t1000\t+\t1\t0\tGT-AG\n")
            second.write_text("chr1\t10\t20\tj1\t1000\t+\t1\t0\tGC-AG\n")
            with self.assertRaisesRegex(ValueError, "conflicting splice motifs"):
                AGGREGATOR.aggregate([first, second], output)

    def test_build_plan_expands_memberships_and_preserves_source_paths(self) -> None:
        fixture_dir = ROOT / "test"
        tracks = VALIDATOR.read_tsv(fixture_dir / "track_manifest.tsv")
        memberships = VALIDATOR.read_tsv(fixture_dir / "aggregation_manifest.tsv")
        plan, metadata, assemblies = VALIDATOR.build_plan(
            tracks,
            memberships,
            fixture_dir,
            {"coverage", "splice_junction"},
            False,
        )
        self.assertEqual(len(plan), 9)
        self.assertEqual(len(metadata), 2)
        self.assertEqual(assemblies, ["GCA_000000001.1"])
        coverage_row = next(row for row in plan if row["track_type"] == "coverage_forward")
        self.assertTrue(Path(coverage_row["track_path"]).is_absolute())
        self.assertEqual(coverage_row["source_path"], "source/run_1.bam")
        self.assertEqual(json.loads(unquote(coverage_row["metadata_json"])), {"tissue": "brain"})

    def test_build_plan_rejects_incompatible_coverage_normalization(self) -> None:
        fixture_dir = ROOT / "test"
        tracks = VALIDATOR.read_tsv(fixture_dir / "track_manifest.tsv")
        memberships = VALIDATOR.read_tsv(fixture_dir / "aggregation_manifest.tsv")
        tracks[1] = dict(tracks[1])
        tracks[1]["normalization_parameters"] = json.dumps(
            {"tool": "deepTools bamCoverage", "bin_size": 50, "normalization": "CPM", "strand": "reverse"}
        )
        with self.assertRaisesRegex(ValueError, "normalization differs"):
            VALIDATOR.build_plan(tracks, memberships, fixture_dir, {"coverage"}, False)

    def test_build_plan_honours_selected_track_type(self) -> None:
        fixture_dir = ROOT / "test"
        tracks = VALIDATOR.read_tsv(fixture_dir / "track_manifest.tsv")
        memberships = VALIDATOR.read_tsv(fixture_dir / "aggregation_manifest.tsv")
        plan, _metadata, _assemblies = VALIDATOR.build_plan(
            tracks,
            memberships,
            fixture_dir,
            {"splice_junction"},
            False,
        )
        self.assertEqual({row["track_type"] for row in plan}, {"splice_junction"})

    def test_build_plan_rejects_aggregate_id_sanitization_collision(self) -> None:
        fixture_dir = ROOT / "test"
        tracks = VALIDATOR.read_tsv(fixture_dir / "track_manifest.tsv")
        memberships = VALIDATOR.read_tsv(fixture_dir / "aggregation_manifest.tsv")
        memberships = [dict(row) for row in memberships]
        memberships[0]["aggregate_id"] = "sample/a"
        memberships[1]["aggregate_id"] = "sample_a"
        with self.assertRaisesRegex(ValueError, "collide after sanitization"):
            VALIDATOR.build_plan(tracks, memberships, fixture_dir, {"coverage"}, False)


if __name__ == "__main__":
    unittest.main()
