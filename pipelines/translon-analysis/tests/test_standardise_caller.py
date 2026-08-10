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
        self.run_standardiser("rpbp", "calls.bed", "chr1\t10\t40\tORF1\t+\t12\n")

    def test_table_callers(self):
        content = "chrom\ttranscript_id\tstart\tend\tframe\tscore\nchr1\tTX1\t10\t40\t0\t5\n"
        for tool in ("ribocode", "orfquant"):
            self.run_standardiser(tool, "calls.tsv", content)

    def test_native_ribotricer_table(self):
        content = "\t".join([
            "ORF_ID", "ORF_type", "status", "phase_score", "read_count",
            "length", "valid_codons", "valid_codons_ratio", "read_density",
            "transcript_id", "transcript_type", "gene_id", "gene_name",
            "gene_type", "chrom", "strand", "start_codon", "profile",
        ]) + "\n"
        content += "\t".join([
            "TX1_101_190_90", "uORF", "translating", "0.91", "25", "90",
            "30", "1.0", "0.28", "TX1", "protein_coding", "G1", "GENE1",
            "protein_coding", "chr1", "+", "ATG", "0,1,2",
        ]) + "\n"
        self.run_standardiser("ribotricer", "ribotricer_translating_ORFs.tsv", content)

    def test_native_ribotish_table(self):
        content = "\t".join([
            "Gid", "Tid", "Symbol", "GeneType", "GenomePos", "StartCodon",
            "Start", "Stop", "TisType", "TISGroup", "TISCounts", "TISPvalue",
            "RiboPvalue", "RiboPStatus", "FisherPvalue", "TISQvalue",
            "FrameQvalue", "FisherQvalue", "AALen", "Seq", "AASeq", "Blocks",
        ]) + "\n"
        content += "\t".join([
            "G1", "TX1", "GENE1", "protein_coding", "chr1:101-220:+", "ATG",
            "0", "120", "Novel", "0", "10", "0.01", "0.001", "T", "0.01",
            "0.02", "0.03", "0.04", "40", "ATG", "M", "101-150,181-220",
        ]) + "\n"
        self.run_standardiser("ribotish", "ribotish.txt", content)

    def test_native_price_table(self):
        content = "Id\tLocation\tType\tp value\tGene\n"
        content += "ENST0001_iORF_1\tchr1+:101-130|151-180\tiORF\t0.001\tGENE1\n"
        self.run_standardiser("price", "price.orfs.tsv", content)

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
