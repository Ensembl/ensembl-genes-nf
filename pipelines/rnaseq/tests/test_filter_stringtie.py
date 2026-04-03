"""
Tests for filter_stringtie.py
"""

import gzip
import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from filter_stringtie import (
    Transcript,
    parse_stringtie_gtf,
    filter_transcripts,
    write_gff3,
    _parse_gtf_attrs,
    _float_attr,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GTF_BASIC = textwrap.dedent("""\
    # StringTie version 2.2.3
    chr1\tStringTie\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "8.500000"; FPKM "1.234"; TPM "2.345";
    chr1\tStringTie\texon\t1000\t2000\t.\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1";
    chr1\tStringTie\texon\t3000\t5000\t.\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1";
    chr1\tStringTie\ttranscript\t7000\t9000\t.\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; cov "1.200000"; FPKM "0.5"; TPM "0.8";
    chr1\tStringTie\texon\t7000\t9000\t.\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1";
    chr2\tStringTie\ttranscript\t500\t800\t.\t+\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; cov "0.5"; FPKM "0.1"; TPM "0.2";
    chr2\tStringTie\texon\t500\t800\t.\t+\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1";
""")

# Single-exon short transcript (should fail min_length and min_exons when set)
GTF_SHORT = textwrap.dedent("""\
    chr1\tStringTie\ttranscript\t1000\t1100\t.\t+\t.\tgene_id "STRG.4"; transcript_id "STRG.4.1"; cov "5.0"; FPKM "1.0"; TPM "1.5";
    chr1\tStringTie\texon\t1000\t1100\t.\t+\t.\tgene_id "STRG.4"; transcript_id "STRG.4.1";
""")


def _write_tmp(content: str, suffix='.gtf') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# _parse_gtf_attrs tests
# ---------------------------------------------------------------------------

class TestParseGtfAttrs:
    def test_basic(self):
        attrs = _parse_gtf_attrs('gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "8.5";')
        assert attrs['gene_id'] == 'STRG.1'
        assert attrs['transcript_id'] == 'STRG.1.1'
        assert attrs['cov'] == '8.5'

    def test_empty(self):
        assert _parse_gtf_attrs('') == {}

    def test_no_trailing_semicolon(self):
        attrs = _parse_gtf_attrs('gene_id "G1"')
        assert attrs['gene_id'] == 'G1'


class TestFloatAttr:
    def test_present(self):
        assert _float_attr({'cov': '8.5'}, 'cov') == 8.5

    def test_missing_returns_default(self):
        assert _float_attr({}, 'cov', 0.0) == 0.0

    def test_invalid_returns_default(self):
        assert _float_attr({'cov': 'N/A'}, 'cov', -1.0) == -1.0


# ---------------------------------------------------------------------------
# parse_stringtie_gtf tests
# ---------------------------------------------------------------------------

class TestParseStringtieGtf:
    def test_parses_transcripts(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        assert len(txs) == 3
        os.unlink(path)

    def test_transcript_attributes(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        tx = next(t for t in txs if t.transcript_id == 'STRG.1.1')
        assert tx.seqname == 'chr1'
        assert tx.start == 1000
        assert tx.end == 5000
        assert tx.strand == '+'
        assert tx.coverage == pytest.approx(8.5)
        assert tx.fpkm == pytest.approx(1.234)
        assert tx.gene_id == 'STRG.1'
        os.unlink(path)

    def test_exons_attached(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        tx = next(t for t in txs if t.transcript_id == 'STRG.1.1')
        assert len(tx.exons) == 2
        assert tx.exons[0] == (1000, 2000)
        assert tx.exons[1] == (3000, 5000)
        os.unlink(path)

    def test_comment_lines_skipped(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        assert all(not tx.seqname.startswith('#') for tx in txs)
        os.unlink(path)

    def test_gz_support(self):
        path = _write_tmp('', suffix='.gtf.gz')
        with gzip.open(path, 'wt') as fh:
            fh.write(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        assert len(txs) == 3
        os.unlink(path)

    def test_minus_strand(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        tx = next(t for t in txs if t.transcript_id == 'STRG.2.1')
        assert tx.strand == '-'
        os.unlink(path)

    def test_length_computed_from_exons(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        tx = next(t for t in txs if t.transcript_id == 'STRG.1.1')
        # exon1 = 2000-1000+1 = 1001, exon2 = 5000-3000+1 = 2001 → 3002
        assert tx.length == 3002
        os.unlink(path)

    def test_n_exons(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        tx = next(t for t in txs if t.transcript_id == 'STRG.1.1')
        assert tx.n_exons == 2
        os.unlink(path)

    def test_empty_file_returns_empty(self):
        path = _write_tmp('# comment\n')
        txs = parse_stringtie_gtf(path)
        assert txs == []
        os.unlink(path)


# ---------------------------------------------------------------------------
# filter_transcripts tests
# ---------------------------------------------------------------------------

class TestFilterTranscripts:
    def _txs(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        os.unlink(path)
        return txs

    def test_coverage_filter(self):
        txs = self._txs()
        # STRG.1.1 cov=8.5, STRG.2.1 cov=1.2, STRG.3.1 cov=0.5
        result = filter_transcripts(txs, min_coverage=2.0, min_length=0, min_exons=0)
        assert len(result) == 1
        assert result[0].transcript_id == 'STRG.1.1'

    def test_length_filter(self):
        path = _write_tmp(GTF_SHORT)
        txs = parse_stringtie_gtf(path)
        os.unlink(path)
        # Short tx: length = 1100-1000+1 = 101 bp
        result = filter_transcripts(txs, min_coverage=0, min_length=200, min_exons=0)
        assert len(result) == 0

    def test_no_filter_keeps_all(self):
        txs = self._txs()
        result = filter_transcripts(txs, min_coverage=0, min_length=0, min_exons=0)
        assert len(result) == 3

    def test_min_exons_filter(self):
        # All STRG.2.1 and STRG.3.1 have 1 exon; STRG.1.1 has 2
        txs = self._txs()
        result = filter_transcripts(txs, min_coverage=0, min_length=0, min_exons=2)
        assert len(result) == 1
        assert result[0].transcript_id == 'STRG.1.1'


# ---------------------------------------------------------------------------
# write_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def _filtered_txs(self):
        path = _write_tmp(GTF_BASIC)
        txs = parse_stringtie_gtf(path)
        os.unlink(path)
        return filter_transcripts(txs, min_coverage=0, min_length=0, min_exons=0)

    def test_writes_header(self):
        txs = self._filtered_txs()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(txs, out.name, 'SAMP1', 'rnaseq_tissue')
        with open(out.name) as fh:
            content = fh.read()
        assert '##gff-version 3' in content
        os.unlink(out.name)

    def test_biotype_set(self):
        txs = self._filtered_txs()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(txs, out.name, 'SAMP1', 'rnaseq_tissue')
        with open(out.name) as fh:
            content = fh.read()
        assert 'biotype=rnaseq_tissue' in content
        os.unlink(out.name)

    def test_gene_transcript_exon_hierarchy(self):
        txs = self._filtered_txs()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(txs, out.name, 'SAMP1', 'rnaseq_tissue')
        with open(out.name) as fh:
            lines = [l for l in fh if '\t' in l]
        features = {l.split('\t')[2] for l in lines}
        assert 'gene' in features
        assert 'transcript' in features
        assert 'exon' in features
        os.unlink(out.name)

    def test_returns_transcript_count(self):
        txs = self._filtered_txs()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        n = write_gff3(txs, out.name, 'SAMP1', 'rnaseq_tissue')
        assert n == len(txs)
        os.unlink(out.name)

    def test_sample_id_in_ids(self):
        txs = self._filtered_txs()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_gff3(txs, out.name, 'MYSAMPLE', 'rnaseq_tissue')
        with open(out.name) as fh:
            content = fh.read()
        assert 'MYSAMPLE' in content
        os.unlink(out.name)
