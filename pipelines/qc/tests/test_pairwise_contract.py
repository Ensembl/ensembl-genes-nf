import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "normalise_pairwise_tsv.py"


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


if __name__ == "__main__":
    unittest.main()
