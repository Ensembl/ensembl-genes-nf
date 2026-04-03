"""
Tests for select_best_targeted.py
"""

import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from select_best_targeted import (
    Gene,
    parse_gff3,
    cluster,
    select_from_cluster,
    write_best_gff3,
    _replace_attr,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _gene(seqname='chr1', start=1000, end=5000, strand='+',
          coverage=90.0, pid=95.0, biotype='cdna_alignment', query_id='NM_001',
          raw_lines=None):
    if raw_lines is None:
        raw_lines = [
            f'{seqname}\texonerate\tgene\t{start}\t{end}\t500\t{strand}\t.'
            f'\tID=gene1;Name={query_id};biotype={biotype};'
            f'coverage={coverage};pid={pid}'
        ]
    return Gene(seqname=seqname, start=start, end=end, strand=strand,
                score=500.0, coverage=coverage, pid=pid,
                query_id=query_id, biotype=biotype, raw_lines=raw_lines)


GFF3_CDNA = textwrap.dedent("""\
    ##gff-version 3
    chr1\texonerate\tgene\t1000\t5000\t500\t+\t.\tID=bt_gene_1;Name=NM_001;biotype=cdna_alignment;coverage=95.0;pid=97.0
    chr1\texonerate\ttranscript\t1000\t5000\t500\t+\t.\tID=bt_tx_1;Parent=bt_gene_1
    chr1\texonerate\texon\t1000\t2000\t.\t+\t.\tID=bt_exon_1;Parent=bt_tx_1
    chr1\texonerate\tgene\t8000\t12000\t400\t-\t.\tID=bt_gene_2;Name=NM_002;biotype=cdna_alignment;coverage=85.0;pid=90.0
    chr1\texonerate\ttranscript\t8000\t12000\t400\t-\t.\tID=bt_tx_2;Parent=bt_gene_2
""")

GFF3_PROTEIN = textwrap.dedent("""\
    ##gff-version 3
    chr1\texonerate\tgene\t1200\t4800\t300\t+\t.\tID=bt_gene_3;Name=P12345;biotype=protein_alignment;coverage=80.0;pid=85.0
    chr1\texonerate\ttranscript\t1200\t4800\t300\t+\t.\tID=bt_tx_3;Parent=bt_gene_3
""")


def _write_gff3(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# parse_gff3 tests
# ---------------------------------------------------------------------------

class TestParseGff3:
    def test_parses_genes(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        assert len(genes) == 2
        os.unlink(path)

    def test_gene_attributes(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        g = genes[0]
        assert g.seqname == 'chr1'
        assert g.start == 1000
        assert g.end == 5000
        assert g.strand == '+'
        assert g.coverage == 95.0
        assert g.pid == 97.0
        assert g.query_id == 'NM_001'
        assert g.biotype == 'cdna_alignment'
        os.unlink(path)

    def test_gene_raw_lines_includes_transcript(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        g = genes[0]
        # raw_lines should include gene + transcript + exon
        assert len(g.raw_lines) >= 2
        os.unlink(path)

    def test_none_returns_empty(self):
        assert parse_gff3(None) == []

    def test_empty_file_returns_empty(self):
        path = _write_gff3('##gff-version 3\n')
        assert parse_gff3(path) == []
        os.unlink(path)


# ---------------------------------------------------------------------------
# cluster tests
# ---------------------------------------------------------------------------

class TestCluster:
    def test_empty_input(self):
        assert cluster([]) == []

    def test_non_overlapping_genes_in_separate_clusters(self):
        g1 = _gene(start=1000, end=2000)
        g2 = _gene(start=5000, end=6000)
        clusters = cluster([g1, g2])
        assert len(clusters) == 2

    def test_overlapping_genes_in_same_cluster(self):
        g1 = _gene(start=1000, end=3000)
        g2 = _gene(start=2500, end=5000)
        clusters = cluster([g1, g2])
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_different_strands_separate(self):
        g1 = _gene(start=1000, end=5000, strand='+')
        g2 = _gene(start=1000, end=5000, strand='-')
        clusters = cluster([g1, g2])
        assert len(clusters) == 2

    def test_different_seqnames_separate(self):
        g1 = _gene(seqname='chr1', start=1000, end=5000)
        g2 = _gene(seqname='chr2', start=1000, end=5000)
        clusters = cluster([g1, g2])
        assert len(clusters) == 2

    def test_touching_but_not_overlapping(self):
        g1 = _gene(start=1000, end=2000)
        g2 = _gene(start=2001, end=4000)
        clusters = cluster([g1, g2])
        assert len(clusters) == 2

    def test_chain_clustering(self):
        # g1-g2 overlap, g2-g3 overlap but g1-g3 don't directly
        g1 = _gene(start=1000, end=2000)
        g2 = _gene(start=1900, end=3000)
        g3 = _gene(start=2900, end=4000)
        clusters = cluster([g1, g2, g3])
        assert len(clusters) == 1


# ---------------------------------------------------------------------------
# select_from_cluster tests
# ---------------------------------------------------------------------------

class TestSelectFromCluster:
    def test_single_analysis_keeps_all(self):
        genes = [_gene(start=1000, end=5000, coverage=90.0, pid=95.0),
                 _gene(start=1500, end=4500, coverage=85.0, pid=90.0)]
        result = select_from_cluster(genes)
        assert len(result) == 2

    def test_mixed_prefers_cdna_over_protein(self):
        cdna    = _gene(biotype='cdna_alignment',    coverage=80.0, pid=80.0)
        protein = _gene(biotype='protein_alignment', coverage=90.0, pid=90.0)
        result = select_from_cluster([cdna, protein])
        assert len(result) == 1
        assert result[0].biotype == 'cdna_alignment'

    def test_mixed_keeps_best_quality_within_analysis(self):
        cdna_good = _gene(biotype='cdna_alignment', coverage=95.0, pid=97.0)
        cdna_bad  = _gene(biotype='cdna_alignment', coverage=60.0, pid=65.0)
        protein   = _gene(biotype='protein_alignment', coverage=85.0, pid=90.0)
        result = select_from_cluster([cdna_good, cdna_bad, protein])
        assert len(result) == 1
        assert result[0].coverage == 95.0

    def test_single_gene_cluster_kept(self):
        g = _gene()
        result = select_from_cluster([g])
        assert result == [g]

    def test_quality_is_coverage_plus_pid(self):
        g1 = _gene(coverage=80.0, pid=80.0)  # quality=160
        g2 = _gene(coverage=90.0, pid=70.0)  # quality=160
        g3 = _gene(coverage=90.0, pid=85.0)  # quality=175
        result = select_from_cluster([g1, g2, g3])
        # Single analysis, so all kept
        assert len(result) == 3


# ---------------------------------------------------------------------------
# write_best_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteBestGff3:
    def test_replaces_biotype_with_best_targeted(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        os.unlink(path)

        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_best_gff3(genes, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert 'biotype=best_targeted' in content
        assert 'biotype=cdna_alignment' not in content
        os.unlink(out.name)

    def test_returns_count(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        n = write_best_gff3(genes, out.name)
        assert n == len(genes)
        os.unlink(out.name)

    def test_preserves_non_gene_lines(self):
        path = _write_gff3(GFF3_CDNA)
        genes = parse_gff3(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_best_gff3(genes, out.name)
        with open(out.name) as fh:
            lines = [l for l in fh if '\t' in l]
        features = {l.split('\t')[2] for l in lines}
        assert 'transcript' in features
        os.unlink(out.name)


# ---------------------------------------------------------------------------
# _replace_attr tests
# ---------------------------------------------------------------------------

class TestReplaceAttr:
    def test_replaces_existing(self):
        result = _replace_attr('ID=foo;biotype=cdna_alignment;pid=95', 'biotype', 'best_targeted')
        assert 'biotype=best_targeted' in result
        assert 'biotype=cdna_alignment' not in result

    def test_adds_if_missing(self):
        result = _replace_attr('ID=foo;pid=95', 'biotype', 'best_targeted')
        assert 'biotype=best_targeted' in result

    def test_other_attrs_preserved(self):
        result = _replace_attr('ID=foo;biotype=cdna_alignment;pid=95', 'biotype', 'best_targeted')
        assert 'ID=foo' in result
        assert 'pid=95' in result
