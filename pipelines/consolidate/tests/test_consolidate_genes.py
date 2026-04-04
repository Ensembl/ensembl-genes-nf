"""
Tests for consolidate_genes.py
"""

import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from consolidate_genes import (
    Transcript,
    parse_gff3,
    cluster,
    select_from_cluster,
    write_consolidated_gff3,
    _replace_attr,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFF3_RNASEQ = textwrap.dedent("""\
    ##gff-version 3
    chr1\trnaseq\tgene\t1000\t5000\t.\t+\t.\tID=rna_gene1;biotype=rnaseq_merged
    chr1\trnaseq\ttranscript\t1000\t5000\t.\t+\t.\tID=rna_tx1;Parent=rna_gene1;biotype=rnaseq_merged
    chr1\trnaseq\texon\t1000\t2000\t.\t+\t.\tID=rna_tx1_exon1;Parent=rna_tx1
    chr1\trnaseq\texon\t3000\t5000\t.\t+\t.\tID=rna_tx1_exon2;Parent=rna_tx1
    chr1\trnaseq\tgene\t8000\t10000\t.\t-\t.\tID=rna_gene2;biotype=rnaseq_merged
    chr1\trnaseq\ttranscript\t8000\t10000\t.\t-\t.\tID=rna_tx2;Parent=rna_gene2;biotype=rnaseq_merged
    chr1\trnaseq\texon\t8000\t10000\t.\t-\t.\tID=rna_tx2_exon1;Parent=rna_tx2
""")

GFF3_AB_INITIO = textwrap.dedent("""\
    ##gff-version 3
    chr1\taugustus\tgene\t1200\t4800\t.\t+\t.\tID=ab_gene1;biotype=ab_initio
    chr1\taugustus\ttranscript\t1200\t4800\t.\t+\t.\tID=ab_tx1;Parent=ab_gene1;biotype=ab_initio
    chr1\taugustus\texon\t1200\t2000\t.\t+\t.\tID=ab_tx1_exon1;Parent=ab_tx1
    chr1\taugustus\texon\t3000\t4800\t.\t+\t.\tID=ab_tx1_exon2;Parent=ab_tx1
    chr2\taugustus\tgene\t100\t500\t.\t+\t.\tID=ab_gene2;biotype=ab_initio
    chr2\taugustus\ttranscript\t100\t500\t.\t+\t.\tID=ab_tx2;Parent=ab_gene2;biotype=ab_initio
    chr2\taugustus\texon\t100\t500\t.\t+\t.\tID=ab_tx2_exon1;Parent=ab_tx2
""")


def _write_tmp(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _tx(seqname='chr1', start=1000, end=5000, strand='+', priority=0,
        tx_id='tx1', gene_id='gene1', biotype='protein_coding'):
    return Transcript(
        seqname=seqname, start=start, end=end, strand=strand,
        tx_id=tx_id, gene_id=gene_id, biotype=biotype,
        priority=priority,
        lines=[
            f'{seqname}\tsrc\ttranscript\t{start}\t{end}\t.\t{strand}\t.'
            f'\tID={tx_id};Parent={gene_id};biotype={biotype}'
        ],
    )


# ---------------------------------------------------------------------------
# _replace_attr tests
# ---------------------------------------------------------------------------

class TestReplaceAttr:
    def test_replaces_existing(self):
        result = _replace_attr('ID=tx1;Parent=old_gene', 'Parent', 'new_gene')
        assert 'Parent=new_gene' in result
        assert 'Parent=old_gene' not in result

    def test_adds_if_missing(self):
        result = _replace_attr('ID=tx1', 'Parent', 'new_gene')
        assert 'Parent=new_gene' in result

    def test_other_attrs_preserved(self):
        result = _replace_attr('ID=tx1;biotype=foo;Parent=old', 'Parent', 'new')
        assert 'ID=tx1' in result
        assert 'biotype=foo' in result


# ---------------------------------------------------------------------------
# parse_gff3 tests
# ---------------------------------------------------------------------------

class TestParseGff3:
    def test_parses_transcripts(self):
        path = _write_tmp(GFF3_RNASEQ)
        txs = parse_gff3(path, priority=2)
        os.unlink(path)
        assert len(txs) == 2

    def test_priority_assigned(self):
        path = _write_tmp(GFF3_RNASEQ)
        txs = parse_gff3(path, priority=2)
        os.unlink(path)
        assert all(t.priority == 2 for t in txs)

    def test_exon_lines_in_tx_lines(self):
        path = _write_tmp(GFF3_RNASEQ)
        txs = parse_gff3(path, priority=2)
        os.unlink(path)
        tx = next(t for t in txs if t.tx_id == 'rna_tx1')
        # lines = [transcript_line, exon1, exon2] = 3
        assert len(tx.lines) == 3

    def test_gene_lines_excluded(self):
        path = _write_tmp(GFF3_RNASEQ)
        txs = parse_gff3(path, priority=2)
        os.unlink(path)
        for tx in txs:
            for line in tx.lines:
                assert '\tgene\t' not in line

    def test_empty_file(self):
        path = _write_tmp('##gff-version 3\n')
        txs = parse_gff3(path, priority=0)
        os.unlink(path)
        assert txs == []

    def test_biotype_parsed(self):
        path = _write_tmp(GFF3_RNASEQ)
        txs = parse_gff3(path, priority=2)
        os.unlink(path)
        assert all(t.biotype == 'rnaseq_merged' for t in txs)


# ---------------------------------------------------------------------------
# cluster tests
# ---------------------------------------------------------------------------

class TestCluster:
    def test_empty(self):
        assert cluster([]) == []

    def test_overlapping_same_strand(self):
        t1 = _tx(start=1000, end=5000)
        t2 = _tx(start=3000, end=7000)
        clusters = cluster([t1, t2])
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_non_overlapping(self):
        t1 = _tx(start=1000, end=2000)
        t2 = _tx(start=5000, end=8000)
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_different_strands_separate(self):
        t1 = _tx(start=1000, end=5000, strand='+', tx_id='a')
        t2 = _tx(start=1000, end=5000, strand='-', tx_id='b')
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_different_seqnames_separate(self):
        t1 = _tx(seqname='chr1', tx_id='a')
        t2 = _tx(seqname='chr2', tx_id='b')
        clusters = cluster([t1, t2])
        assert len(clusters) == 2

    def test_cross_source_cluster(self):
        # rnaseq 1000-5000 and ab_initio 1200-4800 should cluster
        path_rna = _write_tmp(GFF3_RNASEQ)
        path_ab  = _write_tmp(GFF3_AB_INITIO)
        txs = parse_gff3(path_rna, 2) + parse_gff3(path_ab, 5)
        os.unlink(path_rna); os.unlink(path_ab)
        clusters = cluster(txs)
        chr1_plus = [cl for cl in clusters
                     if all(t.seqname == 'chr1' and t.strand == '+' for t in cl)]
        assert len(chr1_plus) == 1
        assert len(chr1_plus[0]) == 2   # one from each source


# ---------------------------------------------------------------------------
# select_from_cluster tests
# ---------------------------------------------------------------------------

class TestSelectFromCluster:
    def test_single_source_keeps_all(self):
        t1 = _tx(start=1000, end=5000, priority=2, tx_id='a')
        t2 = _tx(start=2000, end=6000, priority=2, tx_id='b')
        selected = select_from_cluster([t1, t2])
        assert len(selected) == 2

    def test_higher_priority_wins_in_overlap(self):
        # priority 0 (high quality) overlaps priority 5 (ab initio)
        high = _tx(start=1000, end=5000, priority=0, tx_id='high', gene_id='g_high')
        low  = _tx(start=2000, end=6000, priority=5, tx_id='low',  gene_id='g_low')
        selected = select_from_cluster([high, low])
        ids = [t.tx_id for t in selected]
        assert 'high' in ids
        assert 'low' not in ids

    def test_non_overlapping_lower_priority_kept(self):
        # high priority at 1000-5000; low priority at 8000-10000 (no overlap)
        high = _tx(start=1000, end=5000, priority=0, tx_id='high', gene_id='g1')
        low  = _tx(start=8000, end=10000, priority=5, tx_id='low', gene_id='g2')
        selected = select_from_cluster([high, low])
        ids = [t.tx_id for t in selected]
        assert 'high' in ids
        assert 'low' in ids

    def test_empty_cluster(self):
        assert select_from_cluster([]) == []

    def test_all_same_priority_all_kept(self):
        txs = [_tx(start=1000+i*100, end=2000+i*100, priority=3,
                   tx_id=f'tx{i}', gene_id=f'g{i}') for i in range(3)]
        selected = select_from_cluster(txs)
        assert len(selected) == 3


# ---------------------------------------------------------------------------
# write_consolidated_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteConsolidatedGff3:
    def _selected(self):
        path_rna = _write_tmp(GFF3_RNASEQ)
        path_ab  = _write_tmp(GFF3_AB_INITIO)
        txs = parse_gff3(path_rna, 2) + parse_gff3(path_ab, 5)
        os.unlink(path_rna); os.unlink(path_ab)
        clusters = cluster(txs)
        return [select_from_cluster(cl) for cl in clusters]

    def test_writes_header(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_consolidated_gff3(selected, out.name)
        with open(out.name) as fh:
            assert fh.readline().strip() == '##gff-version 3'
        os.unlink(out.name)

    def test_gene_features_written(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_consolidated_gff3(selected, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\tgene\t' in content
        os.unlink(out.name)

    def test_gene_ids_prefixed(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_consolidated_gff3(selected, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        for line in gene_lines:
            gid = line.split('ID=')[1].split(';')[0]
            assert gid.startswith('consolidated_gene_')
        os.unlink(out.name)

    def test_exon_lines_preserved(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_consolidated_gff3(selected, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\texon\t' in content
        os.unlink(out.name)

    def test_returns_gene_and_tx_counts(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        g_count, tx_count = write_consolidated_gff3(selected, out.name)
        assert g_count > 0
        assert tx_count > 0
        os.unlink(out.name)

    def test_source_gene_in_name_attr(self):
        selected = self._selected()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_consolidated_gff3(selected, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        assert all('Name=' in l for l in gene_lines)
        os.unlink(out.name)
