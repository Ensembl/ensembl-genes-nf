#!/usr/bin/env python3
"""Regression tests for tab-delimited GTF validation."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "validate_gtf.py"
SPEC = importlib.util.spec_from_file_location("validate_gtf", SCRIPT)
assert SPEC and SPEC.loader
validate_gtf_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validate_gtf_module)


class ValidateGtfTests(unittest.TestCase):
    def validate_text(self, text: str) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            gtf = Path(directory) / "input.gtf"
            gtf.write_text(text, encoding="utf-8")
            return validate_gtf_module.validate_gtf(gtf)

    def test_accepts_spaces_in_attributes(self) -> None:
        row = (
            '1\tscallop\texon\t10\t20\t1000\t+\t.\t'
            'gene_id "gene.1"; transcript_id "tx.1"; exon "1";\n'
        )
        self.assertEqual(self.validate_text(row), [])

    def test_ignores_comments_and_blank_lines(self) -> None:
        row = '1\tscallop\texon\t10\t20\t1000\t+\t.\tgene_id "gene.1";\n'
        self.assertEqual(self.validate_text(f"# comment\n\n{row}"), [])

    def test_reports_missing_or_extra_tab_delimited_columns(self) -> None:
        eight_columns = "1\tscallop\texon\t10\t20\t1000\t+\t.\n"
        ten_columns = "1\tscallop\texon\t10\t20\t1000\t+\t.\textra\tfield\n"
        for row, count in ((eight_columns, 8), (ten_columns, 10)):
            with self.subTest(count=count):
                errors = self.validate_text(row)
                self.assertEqual(len(errors), 1)
                self.assertIn(f"found {count}", errors[0])

    def test_rejects_invalid_coordinates_with_line_number(self) -> None:
        row = '1\tscallop\texon\t20\t10\t1000\t+\t.\tgene_id "gene.1";\n'
        errors = self.validate_text(row)
        self.assertEqual(len(errors), 1)
        self.assertIn("line 1", errors[0])
        self.assertIn("1 <= start <= end", errors[0])

    def test_rejects_empty_gtf(self) -> None:
        self.assertEqual(self.validate_text("# comment only\n\n"), ["GTF has no data rows"])


if __name__ == "__main__":
    unittest.main()
