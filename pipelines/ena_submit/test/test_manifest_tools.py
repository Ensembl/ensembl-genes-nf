import csv
import subprocess
import sys
import tempfile
import unittest
from pipelines.ena_submit.bin.poll_webin import extract_accession, receipt_failed
from pathlib import Path


PIPELINE_DIR = Path(__file__).parents[1]
EXPANDER = PIPELINE_DIR / "bin" / "expand_file_manifest.py"
RUN_MANIFEST_BUILDER = PIPELINE_DIR / "bin" / "build_manifest_from_runs.py"


class ExpandManifestTests(unittest.TestCase):
    def test_runs_builder_emits_workflow_compatible_nested_manifest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            input_csv = temp / "runs.csv"
            alignments = temp / "alignments"
            alignments.mkdir()
            (alignments / "ERR1_GCA_1_R1_aln.bam").write_bytes(b"bam")
            input_csv.write_text(
                "run_accession,sample_accession,assembly_accession,taxon_id\n"
                "ERR1,SAMEA1,GCA_1,7955\n"
            )
            outdir = temp / "manifest"
            command = [
                sys.executable, str(RUN_MANIFEST_BUILDER),
                "--input", str(input_csv), "--study", "PRJEB1",
                "--release", "R1", "--mode", "per-run",
                "--file-dir", str(alignments), "--file-ext", "bam",
                "--outdir", str(outdir),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (outdir / "manifest.tsv").open(newline="") as handle:
                manifest = next(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(manifest["files_tsv"].endswith("files/ERR1_GCA_1_R1_aln.tsv"))
            with Path(manifest["files_tsv"]).open(newline="") as handle:
                file_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(file_row["run_accessions"], "ERR1")
            self.assertEqual(file_row["sample_accessions"], "SAMEA1")

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

    def test_poll_receipt_extracts_object_not_submission_accession(self):
        receipt = '<RECEIPT><SUBMISSION accession="ERA1"/><ANALYSIS accession="ERZ1"/></RECEIPT>'
        self.assertEqual(extract_accession(receipt), "ERZ1")
        self.assertTrue(receipt_failed('<RECEIPT success="false"/>'))


if __name__ == "__main__":
    unittest.main()
