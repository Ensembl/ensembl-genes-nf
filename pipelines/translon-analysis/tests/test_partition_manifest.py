import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "make_partition_manifest.py"


class PartitionManifestTest(unittest.TestCase):
    def test_transcriptome_partitions_are_deterministic_and_padded(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "models.bed12"
            bed.write_text("chr2\t100\t200\tTX2\nchr1\t10\t50\tTX1\nchr1\t60\t90\tTX3\n")
            output = Path(tmp) / "manifest.tsv"
            subprocess.run([
                "python3", str(SCRIPT), "--mode", "transcriptome", "--bed12", str(bed),
                "--partitions", "2", "--padding", "5", "--output", str(output)
            ], check=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["contig"] for row in rows], ["chr1", "chr1", "chr2"])
            self.assertEqual({row["partition_id"] for row in rows}, {"0", "1"})
            self.assertEqual(rows[0]["start"], "5")
            self.assertEqual(rows[0]["annotation_load"], "2")
            self.assertEqual(rows[0]["estimated_read_load"], "")

    def test_genome_windows_require_fai_and_include_load_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "empty.bed12"
            fai = Path(tmp) / "reference.fai"
            bed.write_text("")
            fai.write_text("chr1\t25\t0\t25\t26\n")
            output = Path(tmp) / "manifest.tsv"
            subprocess.run([
                "python3", str(SCRIPT), "--mode", "genome", "--bed12", str(bed), "--fai", str(fai),
                "--partitions", "2", "--window-size", "10", "--output", str(output)
            ], check=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(row["start"], row["end"]) for row in rows], [("0", "10"), ("10", "20"), ("20", "25")])
            self.assertIn("estimated_read_load", rows[0])

    def test_transcriptome_gtf_uses_transcript_ids_as_contigs(self):
        with tempfile.TemporaryDirectory() as tmp:
            gtf = Path(tmp) / "transcriptome.gtf"
            gtf.write_text('TX2\ttranslon\tgene\t1\t20\t.\t+\t.\tgene_id "G2"; transcript_id "TX2";\nTX1\ttranslon\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "TX1";\n')
            output = Path(tmp) / "manifest.tsv"
            subprocess.run(["python3", str(SCRIPT), "--mode", "transcriptome", "--gtf", str(gtf), "--output", str(output)], check=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["contig"] for row in rows], ["TX1", "TX2"])


if __name__ == "__main__":
    unittest.main()
