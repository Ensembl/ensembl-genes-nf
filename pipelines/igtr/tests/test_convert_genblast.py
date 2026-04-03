"""Tests for convert_genblast.py"""

from __future__ import annotations

import sys
import os
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from convert_genblast import (
    parse_attrs,
    attrs_to_str,
    load_protein_biotypes,
    parse_genblast_gff,
    passes_filter,
    write_gff3,
    _infer_biotype,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

IMGT_FASTA = textwrap.dedent("""\
    >IG_V_gene|IGHV1-2*02|Homo sapiens
    ATCGATCGATCG
    >TR_J_gene|TRAJ1*01|Homo sapiens
    GCTAGCTAGCTA
    >IG_C_gene|IGHG1*01|Homo sapiens
    TTTTTTTTTTTT
    >UNKNOWNPROT
    AAAAAAAAAA
""")

GENBLAST_GFF = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=IGHV1-2*02-R1-1-A1;Name=IGHV1-2*02;PID=85.50;Coverage=98.30;Note=PID:85.50-Cover:98.30
    chr1\tgenBlastG\tcoding_exon\t1000\t1200\t.\t+\t0\tID=IGHV1-2*02-R1-1-A1-E1;Parent=IGHV1-2*02-R1-1-A1
    chr1\tgenBlastG\tcoding_exon\t1500\t2000\t.\t+\t0\tID=IGHV1-2*02-R1-1-A1-E2;Parent=IGHV1-2*02-R1-1-A1
    chr2\tgenBlastG\ttranscript\t5000\t6000\t62.0\t-\t.\tID=TRAJ1*01-R1-1-A1;Name=TRAJ1*01;PID=62.00;Coverage=70.00;Note=PID:62.00-Cover:70.00
    chr2\tgenBlastG\tcoding_exon\t5000\t6000\t.\t-\t0\tID=TRAJ1*01-R1-1-A1-E1;Parent=TRAJ1*01-R1-1-A1
    chr3\tgenBlastG\ttranscript\t100\t200\t40.0\t+\t.\tID=IGHG1*01-R1-1-A1;Name=IGHG1*01;PID=40.00;Coverage=50.00;Note=PID:40.00-Cover:50.00
    chr3\tgenBlastG\tcoding_exon\t100\t200\t.\t+\t0\tID=IGHG1*01-R1-1-A1-E1;Parent=IGHG1*01-R1-1-A1
""")


@pytest.fixture
def fasta_file(tmp_path):
    f = tmp_path / 'proteins.fa'
    f.write_text(IMGT_FASTA)
    return str(f)


@pytest.fixture
def gff_file(tmp_path):
    f = tmp_path / 'genblast.gff'
    f.write_text(GENBLAST_GFF)
    return str(f)


# ---------------------------------------------------------------------------
# parse_attrs
# ---------------------------------------------------------------------------

class TestParseAttrs:
    def test_simple(self):
        assert parse_attrs('ID=foo;Name=bar') == {'ID': 'foo', 'Name': 'bar'}

    def test_single(self):
        assert parse_attrs('ID=only') == {'ID': 'only'}

    def test_empty(self):
        assert parse_attrs('') == {}

    def test_with_trailing_semicolon(self):
        result = parse_attrs('ID=foo;Name=bar;')
        assert result.get('ID') == 'foo'
        assert result.get('Name') == 'bar'


# ---------------------------------------------------------------------------
# _infer_biotype
# ---------------------------------------------------------------------------

class TestInferBiotype:
    def test_ig_from_raw(self):
        assert _infer_biotype('IG_V_gene', 'X') == 'ig_gene'

    def test_tr_from_raw(self):
        assert _infer_biotype('TR_J_gene', 'X') == 'tr_gene'

    def test_ig_from_name(self):
        assert _infer_biotype('', 'IGHV1-2*02') == 'ig_gene'

    def test_tr_from_name(self):
        assert _infer_biotype('', 'TRAJ1*01') == 'tr_gene'

    def test_unknown_defaults_ig(self):
        assert _infer_biotype('', 'UNKNOWN') == 'ig_gene'


# ---------------------------------------------------------------------------
# load_protein_biotypes
# ---------------------------------------------------------------------------

class TestLoadProteinBiotypes:
    def test_ig_parsed(self, fasta_file):
        bmap = load_protein_biotypes(fasta_file)
        assert bmap['IGHV1-2*02'] == 'ig_gene'

    def test_tr_parsed(self, fasta_file):
        bmap = load_protein_biotypes(fasta_file)
        assert bmap['TRAJ1*01'] == 'tr_gene'

    def test_ig_c_gene_parsed(self, fasta_file):
        bmap = load_protein_biotypes(fasta_file)
        assert bmap['IGHG1*01'] == 'ig_gene'

    def test_header_without_pipe_uses_heuristic(self, fasta_file):
        bmap = load_protein_biotypes(fasta_file)
        # UNKNOWNPROT has no pipe — falls back to name heuristic
        assert 'UNKNOWNPROT' in bmap


# ---------------------------------------------------------------------------
# parse_genblast_gff
# ---------------------------------------------------------------------------

class TestParseGenblastGff:
    def test_transcript_count(self, gff_file):
        txs, exons = parse_genblast_gff(gff_file)
        assert len(txs) == 3

    def test_exon_grouping(self, gff_file):
        txs, exons = parse_genblast_gff(gff_file)
        assert len(exons['IGHV1-2*02-R1-1-A1']) == 2
        assert len(exons['TRAJ1*01-R1-1-A1']) == 1

    def test_pid_attribute_parsed(self, gff_file):
        txs, _ = parse_genblast_gff(gff_file)
        rec = txs['IGHV1-2*02-R1-1-A1']
        assert rec.attrs['PID'] == '85.50'

    def test_coverage_attribute_parsed(self, gff_file):
        txs, _ = parse_genblast_gff(gff_file)
        rec = txs['IGHV1-2*02-R1-1-A1']
        assert rec.attrs['Coverage'] == '98.30'


# ---------------------------------------------------------------------------
# passes_filter
# ---------------------------------------------------------------------------

class TestPassesFilter:
    def setup_method(self):
        from convert_genblast import GffRecord
        self.GffRecord = GffRecord

    def _make_rec(self, pid, cov):
        return self.GffRecord(
            seqname='chr1', source='genBlastG', feature='transcript',
            start=1, end=100, score='.', strand='+', frame='.',
            attrs={'PID': str(pid), 'Coverage': str(cov)},
        )

    def test_passes_both(self):
        assert passes_filter(self._make_rec(85, 95), 70, 80) is True

    def test_fails_pid(self):
        assert passes_filter(self._make_rec(50, 95), 70, 80) is False

    def test_fails_coverage(self):
        assert passes_filter(self._make_rec(85, 60), 70, 80) is False

    def test_exactly_at_threshold(self):
        assert passes_filter(self._make_rec(70, 80), 70, 80) is True


# ---------------------------------------------------------------------------
# write_gff3
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def test_passing_transcripts_written(self, tmp_path, gff_file, fasta_file):
        out = str(tmp_path / 'out.gff3')
        bmap = load_protein_biotypes(fasta_file)
        txs, exons = parse_genblast_gff(gff_file)
        # Default: PID>=70, Coverage>=80 — only IGHV1-2*02 passes
        n = write_gff3(txs, exons, bmap, out, 'test', 70.0, 80.0)
        assert n == 1

    def test_output_has_gene_transcript_exon(self, tmp_path, gff_file, fasta_file):
        out = str(tmp_path / 'out.gff3')
        bmap = load_protein_biotypes(fasta_file)
        txs, exons = parse_genblast_gff(gff_file)
        write_gff3(txs, exons, bmap, out, 'test', 70.0, 80.0)
        content = open(out).read()
        assert '\tgene\t' in content
        assert '\ttranscript\t' in content
        assert '\texon\t' in content

    def test_biotype_in_output(self, tmp_path, gff_file, fasta_file):
        out = str(tmp_path / 'out.gff3')
        bmap = load_protein_biotypes(fasta_file)
        txs, exons = parse_genblast_gff(gff_file)
        write_gff3(txs, exons, bmap, out, 'test', 70.0, 80.0)
        content = open(out).read()
        assert 'biotype=ig_gene' in content

    def test_lower_threshold_includes_tr(self, tmp_path, gff_file, fasta_file):
        out = str(tmp_path / 'out.gff3')
        bmap = load_protein_biotypes(fasta_file)
        txs, exons = parse_genblast_gff(gff_file)
        n = write_gff3(txs, exons, bmap, out, 'test', 60.0, 60.0)
        assert n == 2   # IGHV1-2*02 and TRAJ1*01

    def test_two_exons_written(self, tmp_path, gff_file, fasta_file):
        out = str(tmp_path / 'out.gff3')
        bmap = load_protein_biotypes(fasta_file)
        txs, exons = parse_genblast_gff(gff_file)
        write_gff3(txs, exons, bmap, out, 'test', 70.0, 80.0)
        exon_lines = [l for l in open(out) if '\texon\t' in l]
        assert len(exon_lines) == 2
