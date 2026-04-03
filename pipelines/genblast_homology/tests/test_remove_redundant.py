"""Tests for remove_redundant.py"""

from __future__ import annotations

import sys, os, textwrap
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from remove_redundant import (
    load_genes, cluster_genes, select_from_cluster, write_output, biotype_rank
)


GFF3_MIXED_TIERS = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t.\t+\t.\tID=gene_A;Name=P1;biotype=genblast_3
    chr1\tgenBlastG\ttranscript\t1000\t2000\t.\t+\t.\tID=tx_A;Parent=gene_A;Name=P1;PID=61.0;Coverage=91.0;biotype=genblast_3
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=exon_A;Parent=tx_A
    chr1\tgenBlastG\tgene\t1500\t2500\t.\t+\t.\tID=gene_B;Name=P2;biotype=genblast_1
    chr1\tgenBlastG\ttranscript\t1500\t2500\t.\t+\t.\tID=tx_B;Parent=gene_B;Name=P2;PID=92.0;Coverage=96.0;biotype=genblast_1
    chr1\tgenBlastG\texon\t1500\t2500\t.\t+\t.\tID=exon_B;Parent=tx_B
    chr2\tgenBlastG\tgene\t5000\t6000\t.\t+\t.\tID=gene_C;Name=P3;biotype=genblast_2
    chr2\tgenBlastG\ttranscript\t5000\t6000\t.\t+\t.\tID=tx_C;Parent=gene_C;Name=P3;PID=81.0;Coverage=91.0;biotype=genblast_2
    chr2\tgenBlastG\texon\t5000\t6000\t.\t+\t.\tID=exon_C;Parent=tx_C
""")

GFF3_SAME_TIER_OVERLAPPING = textwrap.dedent("""\
    ##gff-version 3
    chr1\tgenBlastG\tgene\t1000\t2000\t.\t+\t.\tID=gene_A;Name=P1;biotype=genblast_2
    chr1\tgenBlastG\ttranscript\t1000\t2000\t.\t+\t.\tID=tx_A;Parent=gene_A;Name=P1;PID=81.0;Coverage=91.0;biotype=genblast_2
    chr1\tgenBlastG\texon\t1000\t2000\t.\t+\t.\tID=exon_A;Parent=tx_A
    chr1\tgenBlastG\tgene\t1500\t2500\t.\t+\t.\tID=gene_B;Name=P2;biotype=genblast_2
    chr1\tgenBlastG\ttranscript\t1500\t2500\t.\t+\t.\tID=tx_B;Parent=gene_B;Name=P2;PID=82.0;Coverage=92.0;biotype=genblast_2
    chr1\tgenBlastG\texon\t1500\t2500\t.\t+\t.\tID=exon_B;Parent=tx_B
""")


@pytest.fixture
def mixed_file(tmp_path):
    f = tmp_path / 'mixed.gff3'
    f.write_text(GFF3_MIXED_TIERS)
    return str(f)


@pytest.fixture
def same_tier_file(tmp_path):
    f = tmp_path / 'same.gff3'
    f.write_text(GFF3_SAME_TIER_OVERLAPPING)
    return str(f)


class TestBiotypeRank:
    def test_genblast_1_highest(self):
        assert biotype_rank('genblast_1') == 1

    def test_genblast_7_lowest(self):
        assert biotype_rank('genblast_7') == 7

    def test_unknown_gets_99(self):
        assert biotype_rank('unknown') == 99

    def test_ordering(self):
        ranks = [biotype_rank(f'genblast_{i}') for i in range(1, 8)]
        assert ranks == sorted(ranks)


class TestLoadGenes:
    def test_three_genes_loaded(self, mixed_file):
        genes = load_genes(mixed_file)
        assert len(genes) == 3

    def test_biotype_parsed(self, mixed_file):
        genes = load_genes(mixed_file)
        bts = {g.gene_id: g.biotype for g in genes}
        assert bts['gene_A'] == 'genblast_3'
        assert bts['gene_B'] == 'genblast_1'


class TestSelectFromCluster:
    def _gene(self, bt):
        from remove_redundant import Gene
        return Gene('chr1', '+', 1, 100, 'id', 'name', bt, '.', ['line'])

    def test_highest_priority_kept(self):
        cluster = [self._gene('genblast_3'), self._gene('genblast_1')]
        kept = select_from_cluster(cluster)
        assert all(g.biotype == 'genblast_1' for g in kept)

    def test_same_tier_all_kept(self):
        cluster = [self._gene('genblast_2'), self._gene('genblast_2')]
        kept = select_from_cluster(cluster)
        assert len(kept) == 2


class TestWriteOutput:
    def test_high_priority_wins(self, tmp_path, mixed_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(mixed_file)
        clusters = cluster_genes(genes)
        kept, removed = write_output(clusters, out)
        content = open(out).read()
        # gene_B (genblast_1) should be kept; gene_A (genblast_3) removed
        assert 'gene_B' in content
        assert 'gene_A' not in content
        assert removed == 1

    def test_non_overlapping_chr2_kept(self, tmp_path, mixed_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(mixed_file)
        clusters = cluster_genes(genes)
        kept, _ = write_output(clusters, out)
        content = open(out).read()
        assert 'gene_C' in content

    def test_same_tier_both_kept(self, tmp_path, same_tier_file):
        out = str(tmp_path / 'out.gff3')
        genes = load_genes(same_tier_file)
        clusters = cluster_genes(genes)
        kept, removed = write_output(clusters, out)
        assert kept == 2
        assert removed == 0
