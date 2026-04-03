"""
Tests for merge_rnaseq_gff3.py
"""

import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from merge_rnaseq_gff3 import (
    TxRecord,
    parse_gff3_transcripts,
    cluster,
    write_merged_gff3,
    _replace_attr,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFF3_S1 = textwrap.dedent("""\
    ##gff-version 3
    chr1\tStringTie2\tgene\t1000\t5000\t.\t+\t.\tID=S1_rna_gene_STRG.1;biotype=rnaseq_tissue
    chr1\tStringTie2\ttranscript\t1000\t5000\t8.5\t+\t.\tID=S1_rna_tx_STRG.1.1;Parent=S1_rna_gene_STRG.1;biotype=rnaseq_tissue;cov=8.50
    chr1\tStringTie2\texon\t1000\t2000\t.\t+\t.\tID=S1_rna_tx_STRG.1.1_exon_1;Parent=S1_rna_tx_STRG.1.1
    chr1\tStringTie2\texon\t3000\t5000\t.\t+\t.\tID=S1_rna_tx_STRG.1.1_exon_2;Parent=S1_rna_tx_STRG.1.1
    chr1\tStringTie2\tgene\t8000\t10000\t.\t-\t.\tID=S1_rna_gene_STRG.2;biotype=rnaseq_tissue
    chr1\tStringTie2\ttranscript\t8000\t10000\t3.2\t-\t.\tID=S1_rna_tx_STRG.2.1;Parent=S1_rna_gene_STRG.2;biotype=rnaseq_tissue;cov=3.20
    chr1\tStringTie2\texon\t8000\t10000\t.\t-\t.\tID=S1_rna_tx_STRG.2.1_exon_1;Parent=S1_rna_tx_STRG.2.1
""")

GFF3_S2 = textwrap.dedent("""\
    ##gff-version 3
    chr1\tStringTie2\tgene\t1500\t4500\t.\t+\t.\tID=S2_rna_gene_STRG.1;biotype=rnaseq_tissue
    chr1\tStringTie2\ttranscript\t1500\t4500\t6.1\t+\t.\tID=S2_rna_tx_STRG.1.1;Parent=S2_rna_gene_STRG.1;biotype=rnaseq_tissue;cov=6.10
    chr1\tStringTie2\texon\t1500\t2000\t.\t+\t.\tID=S2_rna_tx_STRG.1.1_exon_1;Parent=S2_rna_tx_STRG.1.1
    chr1\tStringTie2\texon\t3000\t4500\t.\t+\t.\tID=S2_rna_tx_STRG.1.1_exon_2;Parent=S2_rna_tx_STRG.1.1
    chr2\tStringTie2\tgene\t100\t500\t.\t+\t.\tID=S2_rna_gene_STRG.2;biotype=rnaseq_tissue
    chr2\tStringTie2\ttranscript\t100\t500\t4.0\t+\t.\tID=S2_rna_tx_STRG.2.1;Parent=S2_rna_gene_STRG.2;biotype=rnaseq_tissue;cov=4.00
    chr2\tStringTie2\texon\t100\t500\t.\t+\t.\tID=S2_rna_tx_STRG.2.1_exon_1;Parent=S2_rna_tx_STRG.2.1
""")


def _write_tmp(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# parse_gff3_transcripts tests
# ---------------------------------------------------------------------------

class TestParseGff3Transcripts:
    def test_parses_transcripts(self):
        path = _write_tmp(GFF3_S1)
        txs = parse_gff3_transcripts(path)
        assert len(txs) == 2
        os.unlink(path)

    def test_transcript_coords(self):
        path = _write_tmp(GFF3_S1)
        txs = parse_gff3_transcripts(path)
        tx = txs[0]
        assert tx.seqname == 'chr1'
        assert tx.start == 1000
        assert tx.end == 5000
        assert tx.strand == '+'
        os.unlink(path)

    def test_exon_lines_attached(self):
        path = _write_tmp(GFF3_S1)
        txs = parse_gff3_transcripts(path)
        tx = txs[0]
        # lines should include transcript + 2 exons
        assert len(tx.lines) == 3
        os.unlink(path)

    def test_gene_lines_excluded(self):
        path = _write_tmp(GFF3_S1)
        txs = parse_gff3_transcripts(path)
        for tx in txs:
            for line in tx.lines:
                assert '\tgene\t' not in line
        os.unlink(path)

    def test_empty_file(self):
        path = _write_tmp('##gff-version 3\n')
        assert parse_gff3_transcripts(path) == []
        os.unlink(path)


# ---------------------------------------------------------------------------
# cluster tests
# ---------------------------------------------------------------------------

class TestCluster:
    def _tx(self, seqname='chr1', start=1000, end=5000, strand='+'):
        return TxRecord(seqname=seqname, start=start, end=end, strand=strand, lines=[])

    def test_empty(self):
        assert cluster([]) == []

    def test_overlapping_same_strand(self):
        t1 = self._tx(start=1000, end=5000)
        t2 = self._tx(start=1500, end=4500)
        clusters = cluster([t1, t2])
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_non_overlapping(self):
        t1 = self._tx(start=1000, end=2000)
        t2 = self._tx(start=5000, end=8000)
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_different_strands_separate(self):
        t1 = self._tx(start=1000, end=5000, strand='+')
        t2 = self._tx(start=1000, end=5000, strand='-')
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_different_seqnames_separate(self):
        t1 = self._tx(seqname='chr1')
        t2 = self._tx(seqname='chr2')
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_cross_sample_clustering(self):
        # One from S1, one from S2, same locus — should be in same cluster
        path1 = _write_tmp(GFF3_S1)
        path2 = _write_tmp(GFF3_S2)
        txs1 = parse_gff3_transcripts(path1)
        txs2 = parse_gff3_transcripts(path2)
        os.unlink(path1)
        os.unlink(path2)
        all_txs = txs1 + txs2
        clusters = cluster(all_txs)
        # chr1+: S1 1000-5000, S2 1500-4500 → same cluster
        chr1_plus = [cl for cl in clusters
                     if all(t.seqname == 'chr1' and t.strand == '+' for t in cl)]
        assert len(chr1_plus) == 1
        assert len(chr1_plus[0]) == 2


# ---------------------------------------------------------------------------
# write_merged_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteMergedGff3:
    def _prepare_clusters(self):
        path1 = _write_tmp(GFF3_S1)
        path2 = _write_tmp(GFF3_S2)
        txs = parse_gff3_transcripts(path1) + parse_gff3_transcripts(path2)
        os.unlink(path1)
        os.unlink(path2)
        return cluster(txs)

    def test_writes_header(self):
        clusters = self._prepare_clusters()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_merged_gff3(clusters, out.name)
        with open(out.name) as fh:
            assert fh.readline().strip() == '##gff-version 3'
        os.unlink(out.name)

    def test_biotype_set_to_merged(self):
        clusters = self._prepare_clusters()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_merged_gff3(clusters, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert 'biotype=rnaseq_merged' in content
        os.unlink(out.name)

    def test_gene_features_written(self):
        clusters = self._prepare_clusters()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        gene_count, tx_count = write_merged_gff3(clusters, out.name)
        assert gene_count > 0
        os.unlink(out.name)

    def test_unified_gene_ids(self):
        clusters = self._prepare_clusters()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_merged_gff3(clusters, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        ids = [l.split('ID=')[1].split(';')[0] for l in gene_lines]
        assert all(id_.startswith('rnaseq_merged_gene_') for id_ in ids)
        os.unlink(out.name)

    def test_exon_lines_preserved(self):
        clusters = self._prepare_clusters()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_merged_gff3(clusters, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\texon\t' in content
        os.unlink(out.name)


# ---------------------------------------------------------------------------
# _replace_attr tests
# ---------------------------------------------------------------------------

class TestReplaceAttr:
    def test_replaces(self):
        result = _replace_attr('ID=old;biotype=rnaseq_tissue', 'biotype', 'rnaseq_merged')
        assert 'biotype=rnaseq_merged' in result
        assert 'biotype=rnaseq_tissue' not in result

    def test_adds_if_missing(self):
        result = _replace_attr('ID=old', 'biotype', 'rnaseq_merged')
        assert 'biotype=rnaseq_merged' in result

    def test_other_attrs_preserved(self):
        result = _replace_attr('ID=old;cov=5.0', 'biotype', 'rnaseq_merged')
        assert 'ID=old' in result
        assert 'cov=5.0' in result
