import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "inflate_sam.awk"


class InflateSamTest(unittest.TestCase):
    def run_inflater(self, sam):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.tsv"
            result = subprocess.run(
                ["awk", "-v", f"manifest={manifest}", "-f", str(SCRIPT)],
                input=sam,
                text=True,
                capture_output=True,
            )
            manifest_text = manifest.read_text() if manifest.exists() else ""
            return result, manifest_text

    def test_expands_read_multiplicity_and_preserves_headers(self):
        sam = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:tx1\tLN:100\n"
            "seq_x1\t0\ttx1\t1\t50\t4M\t*\t0\t0\tACGT\tIIII\n"
            "seq_x2\t0\ttx1\t5\t50\t4M\t*\t0\t0\tTGCA\tIIII\n"
            "seq_x10\t0\ttx1\t9\t50\t4M\t*\t0\t0\tAAAA\tIIII\n"
            "plain\t0\ttx1\t13\t50\t4M\t*\t0\t0\tCCCC\tIIII\n"
        )
        result, manifest = self.run_inflater(sam)
        self.assertEqual(result.returncode, 0, result.stderr)
        records = [line for line in result.stdout.splitlines() if not line.startswith("@")]
        self.assertEqual(len(records), 14)
        self.assertEqual([line.split("\t")[0] for line in records].count("seq_x2"), 2)
        self.assertEqual([line.split("\t")[0] for line in records].count("seq_x10"), 10)
        self.assertIn("inflated_alignments\t14", manifest)

    def test_rejects_malformed_multiplicity_suffix(self):
        sam = "read_xnotanumber\t0\ttx1\t1\t50\t4M\t*\t0\t0\tACGT\tIIII\n"
        result, _ = self.run_inflater(sam)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Malformed multiplicity suffix", result.stderr)


if __name__ == "__main__":
    unittest.main()
