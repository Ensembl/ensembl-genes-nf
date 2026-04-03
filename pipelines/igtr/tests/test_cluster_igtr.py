"""Tests for cluster_igtr.py"""

from __future__ import annotations

import sys
import os
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from cluster_igtr import (
    parse_attrs,
    load_genes,
    cluster_genes,
    write_clustered,
    Gene,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFF3_TWO_OVERLAPPING = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t85.5\t+\t.\tID=gene_A;Name=IGHV1;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=tx_A;Parent=gene_A;Name=IGHV1;PID=85.00;Coverage=95.00;biotype=ig_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=exon_A;Parent=tx_A
    chr1\tgenBlastG\tgene\t1500\t2500\t70.0\t+\t.\tID=gene_B;Name=IGHV2;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1500\t2500\t70.0\t+\t.\tID=tx_B;Parent=gene_B;Name=IGHV2;PID=70.00;Coverage=80.00;biotype=ig_gene
    chr1\tgenBlastG\texon\t1500\t2500\t.\t+\t.\tID=exon_B;Parent=tx_B
""")

GFF3_NON_OVERLAPPING = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t85.5\t+\t.\tID=gene_A;Name=IGHV1;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=tx_A;Parent=gene_A;Name=IGHV1;PID=85.00;Coverage=95.00;biotype=ig_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=exon_A;Parent=tx_A
    chr1\tgenBlastG\tgene\t5000\t6000\t70.0\t+\t.\tID=gene_B;Name=IGHV2;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t5000\t6000\t70.0\t+\t.\tID=tx_B;Parent=gene_B;Name=IGHV2;PID=70.00;Coverage=80.00;biotype=ig_gene
    chr1\tgenBlastG\texon\t5000\t6000\t.\t+\t.\tID=exon_B;Parent=tx_B
""")

GFF3_DIFFERENT_STRANDS = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t85.5\t+\t.\tID=gene_A;Name=IGHV1;biotype=ig_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t85.5\t+\t.\tID=tx_A;Parent=gene_A;Name=IGHV1;PID=85.00;Coverage=95.00;biotype=ig_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=exon_A;Parent=tx_A
    chr1\tgenBlastG\tgene\t1000\t2000\t70.0\t-\t.\tID=gene_B;Name=TRAJ1;biotype=tr_gene
    chr1\tgenBlastG\ttranscript\t1000\t2000\t70.0\t-\t.\tID=tx_B;Parent=gene_B;Name=TRAJ1;PID=70.00;Coverage=80.00;biotype=tr_gene
    chr1\tgenBlastG\texon\t1000\t2000\t.\t-\t.\tID=exon_B;Parent=tx_B
""")


@pytest.fixture
def two_overlap_file(tmp_path):
    f = tmp_path / 'overlap.gff3'
    f.write_text(GFF3_TWO_OVERLAPPING)
    return str(f)


@pytest.fixture
def non_overlap_file(tmp_path):
    f = tmp_path / 'nooverlap.gff3'
    f.write_text(GFF3_NON_OVERLAPPING)
    return str(f)


@pytest.fixture
def diff_strand_file(tmp_path):
    f = tmp_path / 'strands.gff3'
    f.write_text(GFF3_DIFFERENT_STRANDS)
    return str(f)


# ---------------------------------------------------------------------------
# load_genes
# ---------------------------------------------------------------------------

class TestLoadGenes:
    def test_loads_two_genes(self, two_overlap_file):
        genes = load_genes(two_overlap_file)
        assert len(genes) == 2

    def test_pid_parsed(self, two_overlap_file):
        genes = load_genes(two_overlap_file)
        pids = {g.gene_id: g.pid for g in genes}
        assert pids['gene_A'] == pytest.approx(85.0)
        assert pids['gene_B'] == pytest.approx(70.0)

    def test_coverage_parsed(self, two_overlap_file):
        genes = load_genes(two_overlap_file)
        covs = {g.gene_id: g.coverage for g in genes}
        assert covs['gene_A'] == pytest.approx(95.0)

    def test_exon_lines_captured(self, two_overlap_file):
        genes = load_genes(two_overlap_file)
        # Each gene has 3 lines: gene + transcript + exon
        for g in genes:
            assert len(g.lines) == 3


# ---------------------------------------------------------------------------
# Gene.overlaps
# ---------------------------------------------------------------------------

class TestGeneOverlaps:
    def _gene(self, seqname, strand, start, end):
        return Gene(seqname=seqname, strand=strand, start=start, end=end,
                    gene_id='x', name='x', biotype='ig_gene', score='.')

    def test_overlapping(self):
        a = self._gene('chr1', '+', 1000, 2000)
        b = self._gene('chr1', '+', 1500, 2500)
        assert a.overlaps(b)

    def test_adjacent_not_overlapping(self):
        a = self._gene('chr1', '+', 1000, 2000)
        b = self._gene('chr1', '+', 2001, 3000)
        assert not a.overlaps(b)

    def test_different_strand_not_overlapping(self):
        a = self._gene('chr1', '+', 1000, 2000)
        b = self._gene('chr1', '-', 1000, 2000)
        assert not a.overlaps(b)

    def test_different_chrom_not_overlapping(self):
        a = self._gene('chr1', '+', 1000, 2000)
        b = self._gene('chr2', '+', 1000, 2000)
        assert not a.overlaps(b)


# ---------------------------------------------------------------------------
# cluster_genes
# ---------------------------------------------------------------------------

class TestClusterGenes:
    def test_two_overlapping_form_one_cluster(self, two_overlap_file):
        genes = load_genes(two_overlap_file)
        clusters = cluster_genes(genes)
        assert len(clusters) == 1
        assert len(clusters[0]) == 2

    def test_two_non_overlapping_form_two_clusters(self, non_overlap_file):
        genes = load_genes(non_overlap_file)
        clusters = cluster_genes(genes)
        assert len(clusters) == 2

    def test_different_strands_form_two_clusters(self, diff_strand_file):
        genes = load_genes(diff_strand_file)
        clusters = cluster_genes(genes)
        assert len(clusters) == 2


# ---------------------------------------------------------------------------
# write_clustered
# ---------------------------------------------------------------------------

class TestWriteClustered:
    def test_best_gene_selected_by_combined_score(self, tmp_path, two_overlap_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(two_overlap_file)
        clusters = cluster_genes(genes)
        kept, removed = write_clustered(clusters, out)
        assert kept == 1
        assert removed == 1
        # gene_A has PID=85 + Coverage=95 = 180; gene_B has 70+80 = 150
        content = open(out).read()
        assert 'gene_A' in content
        assert 'gene_B' not in content

    def test_non_overlapping_both_kept(self, tmp_path, non_overlap_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(non_overlap_file)
        clusters = cluster_genes(genes)
        kept, removed = write_clustered(clusters, out)
        assert kept == 2
        assert removed == 0

    def test_output_has_gff3_header(self, tmp_path, two_overlap_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(two_overlap_file)
        clusters = cluster_genes(genes)
        write_clustered(clusters, out)
        assert open(out).readline().startswith('##gff-version 3')
