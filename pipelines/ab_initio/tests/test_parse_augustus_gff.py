"""
Tests for parse_augustus_gff.py
"""

import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from parse_augustus_gff import (
    AugGene,
    AugTranscript,
    parse_augustus_gff,
    filter_genes,
    write_gff3,
    _parse_attrs,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

AUG_GFF_BASIC = textwrap.dedent("""\
    # This output was generated with AUGUSTUS (version 3.5.0).
    ##gff-version 3
    chr1\tAUGUSTUS\tgene\t1000\t5000\t.\t+\t.\tID=g1
    chr1\tAUGUSTUS\tmRNA\t1000\t5000\t.\t+\t.\tID=g1.t1;Parent=g1
    chr1\tAUGUSTUS\texon\t1000\t2000\t.\t+\t.\tID=g1.t1.exon1;Parent=g1.t1
    chr1\tAUGUSTUS\texon\t3000\t5000\t.\t+\t.\tID=g1.t1.exon2;Parent=g1.t1
    chr1\tAUGUSTUS\tCDS\t1000\t2000\t.\t+\t0\tID=g1.t1.cds1;Parent=g1.t1
    chr1\tAUGUSTUS\tCDS\t3000\t5000\t.\t+\t0\tID=g1.t1.cds2;Parent=g1.t1
    chr2\tAUGUSTUS\tgene\t200\t800\t.\t-\t.\tID=g2
    chr2\tAUGUSTUS\tmRNA\t200\t800\t.\t-\t.\tID=g2.t1;Parent=g2
    chr2\tAUGUSTUS\texon\t200\t800\t.\t-\t.\tID=g2.t1.exon1;Parent=g2.t1
    chr2\tAUGUSTUS\tCDS\t200\t800\t.\t-\t0\tID=g2.t1.cds1;Parent=g2.t1
""")

# Short gene: 50 bp
AUG_GFF_SHORT = textwrap.dedent("""\
    chr1\tAUGUSTUS\tgene\t1000\t1049\t.\t+\t.\tID=g_short
    chr1\tAUGUSTUS\tmRNA\t1000\t1049\t.\t+\t.\tID=g_short.t1;Parent=g_short
    chr1\tAUGUSTUS\texon\t1000\t1049\t.\t+\t.\tID=g_short.t1.exon1;Parent=g_short.t1
    chr1\tAUGUSTUS\tCDS\t1000\t1049\t.\t+\t0\tID=g_short.t1.cds1;Parent=g_short.t1
""")

# Gene with no explicit exon lines (only CDS) — exons should be derived from CDS
AUG_GFF_NO_EXON = textwrap.dedent("""\
    chr1\tAUGUSTUS\tgene\t1000\t5000\t.\t+\t.\tID=g_noexon
    chr1\tAUGUSTUS\tmRNA\t1000\t5000\t.\t+\t.\tID=g_noexon.t1;Parent=g_noexon
    chr1\tAUGUSTUS\tCDS\t1000\t2000\t.\t+\t0\tID=g_noexon.t1.cds1;Parent=g_noexon.t1
    chr1\tAUGUSTUS\tCDS\t3000\t5000\t.\t+\t0\tID=g_noexon.t1.cds2;Parent=g_noexon.t1
""")


def _write_tmp(content: str, suffix='.gff') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# _parse_attrs tests
# ---------------------------------------------------------------------------

class TestParseAttrs:
    def test_gff3_style(self):
        attrs = _parse_attrs('ID=g1.t1;Parent=g1;biotype=ab_initio')
        assert attrs['ID'] == 'g1.t1'
        assert attrs['Parent'] == 'g1'

    def test_gtf_style(self):
        attrs = _parse_attrs('gene_id "STRG.1"; transcript_id "STRG.1.1";')
        assert attrs['gene_id'] == 'STRG.1'
        assert attrs['transcript_id'] == 'STRG.1.1'

    def test_empty(self):
        assert _parse_attrs('') == {}


# ---------------------------------------------------------------------------
# parse_augustus_gff tests
# ---------------------------------------------------------------------------

class TestParseAugustusGff:
    def test_parses_genes(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        assert len(genes) == 2

    def test_gene_coords(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g1 = next(g for g in genes if g.gene_id == 'g1')
        assert g1.seqname == 'chr1'
        assert g1.start == 1000
        assert g1.end == 5000
        assert g1.strand == '+'

    def test_transcript_attached(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g1 = next(g for g in genes if g.gene_id == 'g1')
        assert len(g1.transcripts) == 1
        assert g1.transcripts[0].tx_id == 'g1.t1'

    def test_exons_attached(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g1 = next(g for g in genes if g.gene_id == 'g1')
        tx = g1.transcripts[0]
        assert len(tx.exons) == 2
        assert tx.exons[0] == (1000, 2000)
        assert tx.exons[1] == (3000, 5000)

    def test_cds_attached(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g1 = next(g for g in genes if g.gene_id == 'g1')
        tx = g1.transcripts[0]
        assert len(tx.cdss) == 2

    def test_minus_strand_gene(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g2 = next(g for g in genes if g.gene_id == 'g2')
        assert g2.strand == '-'

    def test_exons_derived_from_cds_when_absent(self):
        path = _write_tmp(AUG_GFF_NO_EXON)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g = genes[0]
        tx = g.transcripts[0]
        # exons should be derived from CDS
        assert len(tx.exons) == 2
        assert tx.exons[0] == (1000, 2000)

    def test_comment_lines_skipped(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        assert all(not g.seqname.startswith('#') for g in genes)

    def test_empty_file(self):
        path = _write_tmp('# comment\n')
        genes = parse_augustus_gff(path)
        os.unlink(path)
        assert genes == []


# ---------------------------------------------------------------------------
# filter_genes tests
# ---------------------------------------------------------------------------

class TestFilterGenes:
    def _genes(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        return genes

    def test_keep_all_above_threshold(self):
        genes = self._genes()
        result = filter_genes(genes, min_gene_length=100)
        # g1: 5000-1000+1=4001, g2: 800-200+1=601 — both > 100
        assert len(result) == 2

    def test_filter_short_gene(self):
        path = _write_tmp(AUG_GFF_SHORT)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        result = filter_genes(genes, min_gene_length=100)
        assert len(result) == 0

    def test_zero_threshold_keeps_all(self):
        genes = self._genes()
        assert len(filter_genes(genes, 0)) == len(genes)

    def test_gene_length_property(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        g1 = next(g for g in genes if g.gene_id == 'g1')
        assert g1.length == 4001   # 5000 - 1000 + 1


# ---------------------------------------------------------------------------
# write_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def _genes(self):
        path = _write_tmp(AUG_GFF_BASIC)
        genes = parse_augustus_gff(path)
        os.unlink(path)
        return filter_genes(genes, min_gene_length=0)

    def test_writes_header(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            assert fh.readline().strip() == '##gff-version 3'
        os.unlink(out.name)

    def test_biotype_ab_initio(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert 'biotype=ab_initio' in content
        os.unlink(out.name)

    def test_gene_ids_prefixed(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        for line in gene_lines:
            gid = line.split('ID=')[1].split(';')[0]
            assert gid.startswith('ab_initio_gene_')
        os.unlink(out.name)

    def test_exon_features_written(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\texon\t' in content
        os.unlink(out.name)

    def test_cds_features_written(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\tCDS\t' in content
        os.unlink(out.name)

    def test_returns_gene_tx_counts(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        g_count, tx_count = write_gff3(genes, out.name)
        assert g_count == 2
        assert tx_count == 2
        os.unlink(out.name)

    def test_source_gene_id_in_name(self):
        genes = self._genes()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(genes, out.name)
        with open(out.name) as fh:
            content = fh.read()
        # Original IDs should appear as Name= attributes
        assert 'Name=g1' in content
        os.unlink(out.name)
