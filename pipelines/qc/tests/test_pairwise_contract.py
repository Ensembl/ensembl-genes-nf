import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "normalise_pairwise_tsv.py"
PAIRWISE_QC = Path(__file__).parents[1] / "bin" / "pairwise_qc.py"
TESTDATA = Path(__file__).parents[1] / "testdata" / "pairwise"


class PairwiseContractTest(unittest.TestCase):
    def test_legacy_headers_are_normalised_to_source_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "legacy.tsv"
            output = root / "normalised.tsv"
            source.write_text(
                "ensembl_gene_id\tcat_gene_id\tens_to_cat_concordance_rate\n"
                "a1\tb1\t1.0\n",
                encoding="utf-8",
            )
            subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output)],
                check=True,
            )
            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines()[0],
                "source_a_gene_id\tsource_b_gene_id\tsource_a_to_source_b_concordance_rate",
            )

    def test_bed12_normalisation_is_structural_and_has_no_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "models.tsv"
            subprocess.run(
                [sys.executable, str(PAIRWISE_QC), "normalise", "--input", str(TESTDATA / "source_a.bed12"),
                 "--source", "source_a", "--pair-id", "bed_pair", "--assembly", "fixture",
                 "--format", "bed12", "--output", str(output)],
                check=True,
            )
            lines = output.read_text(encoding="utf-8").splitlines()
            self.assertIn("source_format", lines[0])
            self.assertIn("\tbed12\tvalid\t", lines[1])
            self.assertTrue(lines[1].endswith("\t\t"))

    def test_bed12_conversion_marks_synthetic_gene_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "converted.gtf"
            subprocess.run(
                [sys.executable, str(PAIRWISE_QC), "bed12gtf", "--input", str(TESTDATA / "source_a.bed12"), "--output", str(output)],
                check=True,
            )
            text = output.read_text(encoding="utf-8")
            self.assertIn('gene_id "synthetic_gene:lr_tx_1"', text)

    def test_manifest_validation_rejects_unknown_format(self):
        result = subprocess.run(
            [sys.executable, str(PAIRWISE_QC), "validate", "--pair-id", "p", "--assembly", "a",
             "--source-a", "a", "--source-b", "b", "--source-a-path", str(TESTDATA / "source_a.gff3"),
             "--source-b-path", str(TESTDATA / "source_b.gff3"), "--source-a-format", "wat",
             "--source-b-format", "gff3", "--output", "/tmp/invalid-pair.tsv"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_manifest_validation_rejects_mismatched_assemblies(self):
        result = subprocess.run(
            [sys.executable, str(PAIRWISE_QC), "validate", "--pair-id", "p", "--assembly", "a",
             "--source-a", "a", "--source-b", "b", "--source-a-path", str(TESTDATA / "source_a.gff3"),
             "--source-b-path", str(TESTDATA / "source_b.gff3"), "--source-a-format", "gff3",
             "--source-b-format", "gff3", "--source-a-assembly", "a", "--source-b-assembly", "b",
             "--output", "/tmp/invalid-assembly.tsv"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
