import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "make_samplesheet.py"


class MakeSamplesheetTest(unittest.TestCase):
    def test_discovers_paired_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "riboseq"
            (root / "sampleA").mkdir(parents=True)
            for name in ("sampleA.Aligned.toTranscriptome.out.bam", "sampleA.Aligned.sortedByCoord.out.bam"):
                path = root / "sampleA" / name
                path.touch()
                Path(f"{path}.bai").touch()
            output = Path(tmp) / "samplesheet.tsv"
            subprocess.run(["python3", str(SCRIPT), "--riboseq-outdir", str(root), "--output", str(output)], check=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sample_id"], "sampleA")
            self.assertTrue(rows[0]["transcriptome_bai"])
            self.assertTrue(rows[0]["genome_bai"])


if __name__ == "__main__":
    unittest.main()
