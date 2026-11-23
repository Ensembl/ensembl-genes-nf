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
ANNOTATION_MANIFEST_BUILDER = PIPELINE_DIR / "bin" / "build_manifest_from_rnaseq_annotation.py"


class ExpandManifestTests(unittest.TestCase):
    def test_annotation_builder_preserves_alternate_alignment_path_spelling(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp = Path(temp_dir)
            assembly = temp / "danio_rerio" / "GCA_1"
            rnaseq = assembly / "rnaseq"
            rnaseq.mkdir(parents=True)
            fasta = assembly / "GCA_1_genomic.fna"
            fasta.write_text(">chr1\nACGT\n")
            report = assembly / "GCA_1_assembly_report.txt"
            report.write_text("# Sequence-Name\tSequence-Role\n")

            alignment_target = temp / "alignment_target"
            (alignment_target / "output").mkdir(parents=True)
            bam = alignment_target / "output" / "ERR1_Aligned.sortedByCoord.out.bam"
            bam.write_bytes(b"bam")
            alignment_alias = temp / "nfs_alignment"
            alignment_alias.symlink_to(alignment_target, target_is_directory=True)

            runs_csv = rnaseq / "danio_rerio.csv"
            runs_csv.write_text("sample_1\tERR1\n")
            outdir = temp / "manifest"
            command = [
                sys.executable, str(ANNOTATION_MANIFEST_BUILDER),
                "--rnaseq-dir", str(rnaseq),
                "--alignment-dir", str(alignment_alias),
                "--runs-csv", str(runs_csv),
                "--reference-fasta", str(fasta),
                "--assembly-report", str(report),
                "--assembly-accession", "GCA_1",
                "--last-geneset-update", "2026-09",
                "--species", "danio_rerio",
                "--study", "PRJEB1",
                "--file-format", "bam",
                "--no-bam-header",
                "--outdir", str(outdir),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (outdir / "files.tsv").open(newline="") as handle:
                file_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                file_row["file_path"],
                str(alignment_alias / "output" / "ERR1_Aligned.sortedByCoord.out.bam"),
            )

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
