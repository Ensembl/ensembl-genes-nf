import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "standardise_caller.py"


class StandardiseCallerTest(unittest.TestCase):
    def run_standardiser(self, tool, raw_name, raw_content, gtf_content=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw"
            raw.mkdir()
            (raw / raw_name).write_text(raw_content)
            gtf = None
            if gtf_content:
                gtf = root / "annotation.gff3"
                gtf.write_text(gtf_content)
            tsv = root / "S1.tsv"
            bed = root / "S1.bed12"
            subprocess.run([
                sys.executable, str(SCRIPT), "--raw", str(raw), "--tool", tool,
                "--sample", "S1", "--output", str(tsv), "--bed12", str(bed),
                *( ["--gtf", str(gtf)] if gtf else [] ),
            ], check=True)
            with tsv.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["tool"], tool)
            self.assertEqual(rows[0]["chrom"], "chr1")
            self.assertEqual(len(bed.read_text().splitlines()[0].split("\t")), 12)

    def test_bed_callers(self):
        self.run_standardiser("ribotaper", "calls.bed", "chr1\t10\t40\tORF1\t+\t12\n")
        self.run_standardiser("rpbp", "calls.bed", "chr1\t10\t40\tORF1\t+\t12\n")

    def test_table_callers(self):
        content = "chrom\ttranscript_id\tstart\tend\tframe\tscore\nchr1\tTX1\t10\t40\t0\t5\n"
        for tool in ("ribocode", "ribotricer", "orfquant"):
            self.run_standardiser(tool, "calls.tsv", content)

    def test_transcript_coordinates_project_to_genome_blocks(self):
        content = "transcript_id\tstart\tend\tframe\tscore\nTX1\t0\t15\t0\t5\n"
        gff = "\n".join([
            "##gff-version 3",
            "chr1\tTEST\tmRNA\t101\t130\t.\t+\t.\tID=transcript:TX1",
            "chr1\tTEST\texon\t101\t110\t.\t+\t.\tParent=transcript:TX1",
            "chr1\tTEST\texon\t121\t130\t.\t+\t.\tParent=transcript:TX1",
        ])
        self.run_standardiser("ribocode", "calls.tsv", content, gff)


if __name__ == "__main__":
    unittest.main()
