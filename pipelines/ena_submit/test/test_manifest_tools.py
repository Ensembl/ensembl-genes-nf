import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PIPELINE_DIR = Path(__file__).parents[1]
EXPANDER = PIPELINE_DIR / "bin" / "expand_file_manifest.py"


class ExpandManifestTests(unittest.TestCase):
    def test_missing_alignment_is_rejected_unless_stub_check_is_explicitly_skipped(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            files_tsv = temp / "files.tsv"
            analysis_tsv = temp / "analysis.tsv"
            output_tsv = temp / "expanded.tsv"
            files_tsv.write_text(
                "file_path\tfile_type\trun_accession\n"
                "/does/not/exist.bam\tbam\tSRR000001\n"
            )
            analysis_tsv.write_text(
                "analysis_alias\ttitle\tassembly_accession\n"
                "example\tExample\tGCA_000001405.28\n"
            )

            command = [
                sys.executable,
                str(EXPANDER),
                "--files-tsv",
                str(files_tsv),
                "--analysis-tsv",
                str(analysis_tsv),
                "--analysis-id",
                "example",
                "--project-alias",
                "prj_example",
                "--assembly",
                "GCA_000001405.28",
                "--release",
                "release",
                "--out",
                str(output_tsv),
            ]
            rejected = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("missing alignment file", rejected.stderr)

            accepted = subprocess.run(
                command + ["--skip-file-check"], capture_output=True, text=True
            )
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            with output_tsv.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["run_accession"], "SRR000001")


if __name__ == "__main__":
    unittest.main()
