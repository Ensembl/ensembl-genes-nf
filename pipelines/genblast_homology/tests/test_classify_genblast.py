"""Tests for classify_genblast.py"""

from __future__ import annotations

import sys, os, textwrap
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from classify_genblast import classify, parse_genblast_gff, write_classified_gff3, TIERS


GENBLAST_GFF = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\ttranscript\t1000\t2000\t95.0\t+\t.\tID=PROT1-R1-1-A1;Name=PROT1;PID=92.00;Coverage=96.00
    chr1\tgenBlastG\tcoding_exon\t1000\t2000\t.\t+\t0\tID=PROT1-R1-1-A1-E1;Parent=PROT1-R1-1-A1
    chr2\tgenBlastG\ttranscript\t3000\t4000\t80.0\t-\t.\tID=PROT2-R1-1-A1;Name=PROT2;PID=65.00;Coverage=91.00
    chr2\tgenBlastG\tcoding_exon\t3000\t4000\t.\t-\t0\tID=PROT2-R1-1-A1-E1;Parent=PROT2-R1-1-A1
    chr3\tgenBlastG\ttranscript\t5000\t6000\t50.0\t+\t.\tID=PROT3-R1-1-A1;Name=PROT3;PID=25.00;Coverage=55.00
    chr3\tgenBlastG\tcoding_exon\t5000\t6000\t.\t+\t0\tID=PROT3-R1-1-A1-E1;Parent=PROT3-R1-1-A1
    chr4\tgenBlastG\ttranscript\t7000\t8000\t20.0\t+\t.\tID=PROT4-R1-1-A1;Name=PROT4;PID=10.00;Coverage=20.00
    chr4\tgenBlastG\tcoding_exon\t7000\t8000\t.\t+\t0\tID=PROT4-R1-1-A1-E1;Parent=PROT4-R1-1-A1
""")


@pytest.fixture
def gff_file(tmp_path):
    f = tmp_path / 'gb.gff'
    f.write_text(GENBLAST_GFF)
    return str(f)


class TestClassify:
    def test_tier1(self):
        assert classify(92, 96) == 'genblast_1'

    def test_tier2(self):
        assert classify(81, 91) == 'genblast_2'

    def test_tier3(self):
        assert classify(61, 91) == 'genblast_3'

    def test_tier4(self):
        assert classify(41, 91) == 'genblast_4'

    def test_tier5(self):
        assert classify(21, 81) == 'genblast_5'

    def test_tier6(self):
        assert classify(21, 61) == 'genblast_6'

    def test_tier7_low_coverage(self):
        assert classify(50, 30) == 'genblast_7'

    def test_tier7_all_zero(self):
        assert classify(0, 0) == 'genblast_7'

    def test_exact_boundary_tier1(self):
        assert classify(90, 95) == 'genblast_1'


class TestParseGenblastGff:
    def test_count(self, gff_file):
        txs, exons = parse_genblast_gff(gff_file)
        assert len(txs) == 4

    def test_pid_parsed(self, gff_file):
        txs, _ = parse_genblast_gff(gff_file)
        assert txs['PROT1-R1-1-A1'].attrs['PID'] == '92.00'


class TestWriteClassifiedGff3:
    def test_all_pass_at_low_threshold(self, tmp_path, gff_file):
        out = str(tmp_path / 'out.gff3')
        txs, exons = parse_genblast_gff(gff_file)
        counts = write_classified_gff3(txs, exons, out, 'test', 0.0, 0.0)
        assert sum(counts.values()) == 4

    def test_filter_removes_below_threshold(self, tmp_path, gff_file):
        out = str(tmp_path / 'out.gff3')
        txs, exons = parse_genblast_gff(gff_file)
        counts = write_classified_gff3(txs, exons, out, 'test', 30.0, 50.0)
        # PROT3 (PID=25) fails min_pid=30; PROT4 (PID=10,Cov=20) fails both
        assert sum(counts.values()) == 2

    def test_tier_distribution(self, tmp_path, gff_file):
        out = str(tmp_path / 'out.gff3')
        txs, exons = parse_genblast_gff(gff_file)
        counts = write_classified_gff3(txs, exons, out, 'test', 0.0, 0.0)
        assert counts.get('genblast_1', 0) == 1  # PROT1: cov=96 pid=92
        assert counts.get('genblast_3', 0) == 1  # PROT2: cov=91 pid=65

    def test_gene_transcript_exon_present(self, tmp_path, gff_file):
        out = str(tmp_path / 'out.gff3')
        txs, exons = parse_genblast_gff(gff_file)
        write_classified_gff3(txs, exons, out, 'test', 0.0, 0.0)
        content = open(out).read()
        assert '\tgene\t' in content
        assert '\ttranscript\t' in content
        assert '\texon\t' in content

    def test_biotype_in_attrs(self, tmp_path, gff_file):
        out = str(tmp_path / 'out.gff3')
        txs, exons = parse_genblast_gff(gff_file)
        write_classified_gff3(txs, exons, out, 'test', 0.0, 0.0)
        content = open(out).read()
        assert 'biotype=genblast_' in content
