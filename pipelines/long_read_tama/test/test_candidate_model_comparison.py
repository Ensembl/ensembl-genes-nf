#!/usr/bin/env python3
"""Unit tests for the backend-neutral canonical model comparison."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "bin" / "compare_candidate_models.py"
RESOURCE_SCRIPT = ROOT / "bin" / "summarise_backend_resources.py"


BED_A = "chr1\t0\t20\ta\t0\t+\t0\t20\t0\t2\t5,5,\t0,15,\n"
BED_B = "chr1\t0\t25\tb\t0\t+\t0\t25\t0\t2\t5,5,\t0,20,\n"


class CandidateModelComparisonTest(unittest.TestCase):
    def test_summary_and_pairwise_structure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left = root / "tama_combined_models.bed"
            right = root / "isoquant_combined_models.bed"
            left.write_text(BED_A)
            right.write_text(BED_B)
            summary = root / "summary.tsv"
            report = root / "summary.json"
            manifest = root / "manifest.tsv"
            subprocess.run([
                "python3", str(SCRIPT), str(summary), str(report), str(left), str(right),
                "--manifest", str(manifest), "--stage", "cohort",
            ], check=True)

            rows = summary.read_text().splitlines()
            self.assertEqual(len(rows), 3)
            self.assertTrue(any(row.startswith("tama\t") for row in rows[1:]))
            self.assertTrue(any(row.startswith("isoquant\t") for row in rows[1:]))
            payload = json.loads(report.read_text())
            self.assertEqual(len(payload["backends"]), 2)
            self.assertEqual(len(payload["pairwise"]), 1)
            self.assertEqual(payload["pairwise"][0]["shared_unique_intron_chains"], 0)
            manifest_rows = manifest.read_text().splitlines()
            self.assertEqual(len(manifest_rows), 3)
            self.assertIn("cohort", manifest_rows[1])
            self.assertIn("STRUCTURAL_ONLY", manifest_rows[1])

    def test_resource_summary_groups_trace_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            trace = root / "trace.tsv"
            trace.write_text(
                "name\tstatus\trealtime\tpeak_rss\n"
                "RUN_ISOQUANT:ISOQUANT\tCOMPLETED\t1m 2s\t4 GB\n"
                "MERGE_LONG_READ_MODELS:TAMA_MERGE\tFAILED\t3s\t512 MB\n"
                "RUN_BAMBU:BAMBU_DISCOVERY\tFAILED\t2s\t1 GB\n"
            )
            summary = root / "resources.tsv"
            report = root / "resources.json"
            subprocess.run(["python3", str(RESOURCE_SCRIPT), str(trace), str(summary), str(report)], check=True)
            text = summary.read_text()
            self.assertIn("isoquant\t1\t1", text)
            self.assertIn("tama\t1\t0\t1", text)
            self.assertIn("bambu\t1\t0\t1", text)


if __name__ == "__main__":
    unittest.main()
