#!/usr/bin/env python3
"""Regression tests for native transcript shard counting."""

from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "audit_native_shard_models.py"


class GtfShardAuditTest(unittest.TestCase):
    def audit(self, root: Path, name: str, contents: str) -> dict[str, str]:
        model = root / f"{name}.gtf"
        mergeable = root / f"{name}.mergeable.gtf"
        counts = root / f"{name}.counts.tsv"
        model.write_text(contents)
        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                str(model),
                "gtf",
                "flair",
                "ACC1",
                name,
                str(mergeable),
                str(counts),
            ],
            check=True,
        )
        with counts.open() as handle:
            return next(csv.DictReader(handle, delimiter="\t"))

    def test_gene_and_exon_rows_do_not_count_as_transcripts(self) -> None:
        gtf = (
            'chr1\ttool\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n'
            'chr1\ttool\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\ttool\texon\t1\t50\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            row = self.audit(Path(directory), "mixed", gtf)
            self.assertEqual(row["model_count"], "1")
            self.assertEqual(row["merge_status"], "MERGEABLE")

    def test_gtf_without_transcript_rows_is_an_empty_shard(self) -> None:
        gtf = (
            'chr1\ttool\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n'
            'chr1\ttool\texon\t1\t50\t.\t+\t.\tgene_id "g1";\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            row = self.audit(Path(directory), "empty", gtf)
            self.assertEqual(row["model_count"], "0")
            self.assertEqual(row["merge_status"], "EMPTY")


if __name__ == "__main__":
    unittest.main()
