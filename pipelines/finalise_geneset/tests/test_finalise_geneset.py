"""
Tests for finalise_geneset bin scripts.

Covers:
  filter_geneset.py   — ORF length, intron size, whole-gene removal
  detect_pseudogenes.py — repeat coverage, frameshifted introns
  detect_readthrough.py — spanning two gene CDSs
  select_canonical.py  — longest CDS, span tie-break, attribute set
  flag_selenoproteins.py — name matching, NO_FILE passthrough
"""

from __future__ import annotations

import os
import sys
import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Path setup: allow importing from the sibling bin/ directory
# ---------------------------------------------------------------------------
BIN = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, BIN)


# ===========================================================================
# ── filter_geneset ──────────────────────────────────────────────────────────
# ===========================================================================

from filter_geneset import (
    GffRecord as FGffRecord,
    _parse_attrs as fg_parse_attrs,
    _set_attr as fg_set_attr,
    cds_length as fg_cds_length,
    min_intron_size_for_tx,
    parse_gff3 as fg_parse_gff3,
    should_remove_transcript,
    write_filtered_gff3,
)


def make_filter_gff3(tmp_path: Path, content: str) -> Path:
    p = tmp_path / 'input.gff3'
    p.write_text(textwrap.dedent(content))
    return p


class TestFilterGff3Parsing:
    def test_parse_returns_gene(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;biotype=protein_coding
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t1
        """)
        _, genes, transcripts, children, gene_order, tx_order = fg_parse_gff3(str(gff))
        assert 'g1' in genes
        assert 't1' in transcripts
        assert gene_order == ['g1']
        assert tx_order['g1'] == ['t1']
        cds_recs = [r for r in children['t1'] if r.feature == 'CDS']
        assert len(cds_recs) == 1

    def test_parse_multiexon(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t1
            chr1\t.\texon\t600\t2000\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t1
            chr1\t.\tCDS\t600\t2000\t.\t+\t0\tParent=t1
        """)
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert fg_cds_length('t1', children) == 500 + 1401


class TestCdsLength:
    def test_single_cds(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t1\t300\t.\t+\t0\tParent=t1
        """)
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert fg_cds_length('t1', children) == 300

    def test_multi_cds(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t100\t200\t.\t+\t0\tParent=t1
            chr1\t.\tCDS\t300\t400\t.\t+\t0\tParent=t1
        """)
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert fg_cds_length('t1', children) == 101 + 101


class TestMinIntronSize:
    def test_single_exon_returns_none(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t1000\t.\t+\t.\tParent=t1
        """)
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert min_intron_size_for_tx('t1', children) is None

    def test_large_intron(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t1
            chr1\t.\texon\t600\t2000\t.\t+\t.\tParent=t1
        """)
        # Intron: 501..599 = 99 bp
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert min_intron_size_for_tx('t1', children) == 99

    def test_tiny_intron(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t1
            chr1\t.\texon\t504\t2000\t.\t+\t.\tParent=t1
        """)
        # Intron: 501..503 = 3 bp
        _, _, _, children, _, _ = fg_parse_gff3(str(gff))
        assert min_intron_size_for_tx('t1', children) == 3


class TestShouldRemoveTranscript:
    """Unit tests for the removal decision logic."""

    def _make_children(self, cds_ranges, exon_ranges):
        """Build a minimal children dict for tx_id='t1'."""
        from filter_geneset import GffRecord as R
        recs = []
        for s, e in exon_ranges:
            recs.append(R('c', '.', 'exon', s, e, '.', '+', '.', 'Parent=t1', ''))
        for s, e in cds_ranges:
            recs.append(R('c', '.', 'CDS', s, e, '.', '+', '0', 'Parent=t1', ''))
        return {'t1': recs}

    def test_short_orf_removed(self):
        # CDS = 60 bp → 20 aa < 100 aa threshold
        children = self._make_children([(1, 60)], [(1, 60)])
        removed, reason = should_remove_transcript('t1', children, 100, 10)
        assert removed
        assert 'short_orf' in reason

    def test_long_orf_kept(self):
        # CDS = 900 bp → 300 aa
        children = self._make_children([(1, 900)], [(1, 900)])
        removed, _ = should_remove_transcript('t1', children, 100, 10)
        assert not removed

    def test_tiny_intron_removed(self):
        # Two exons with 3 bp intron (exon1: 1-500, exon2: 504-1000)
        # CDS is large enough
        children = self._make_children(
            [(1, 500), (504, 1000)],   # CDS spans both exons, large
            [(1, 500), (504, 1000)],
        )
        removed, reason = should_remove_transcript('t1', children, 10, 10)
        assert removed
        assert 'tiny_intron' in reason

    def test_normal_intron_kept(self):
        # Exons: 1-500, 600-1000 → intron size = 99 bp
        # CDS: 1-500 + 600-1000 = 901 bp → 300 aa
        children = self._make_children(
            [(1, 500), (600, 1000)],
            [(1, 500), (600, 1000)],
        )
        removed, _ = should_remove_transcript('t1', children, 100, 10)
        assert not removed


class TestWriteFilteredGff3:
    def test_gene_removed_when_all_transcripts_removed(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t100\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t1\t60\t.\t+\t0\tParent=t1
        """)
        # CDS = 60 bp < 100*3 = 300 bp → removed
        _, genes, transcripts, children, gene_order, tx_order = fg_parse_gff3(str(gff))

        removed_txs = {'t1'}
        removed_genes = {'g1'}
        out = tmp_path / 'out.gff3'
        write_filtered_gff3(
            str(out), ['##gff-version 3\n'],
            genes, transcripts, children, gene_order, tx_order,
            removed_txs, removed_genes,
        )
        content = out.read_text()
        assert 'g1' not in content
        assert 't1' not in content

    def test_kept_gene_written(self, tmp_path):
        gff = make_filter_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t1000\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t1\t900\t.\t+\t0\tParent=t1
        """)
        _, genes, transcripts, children, gene_order, tx_order = fg_parse_gff3(str(gff))
        out = tmp_path / 'out.gff3'
        write_filtered_gff3(
            str(out), ['##gff-version 3\n'],
            genes, transcripts, children, gene_order, tx_order,
            set(), set(),
        )
        content = out.read_text()
        assert 'g1' in content
        assert 't1' in content


# ===========================================================================
# ── detect_pseudogenes ──────────────────────────────────────────────────────
# ===========================================================================

from detect_pseudogenes import (
    GffRecord as PGffRecord,
    all_introns_frameshifted,
    classify_gene,
    covered_bases,
    exon_count,
    has_cds,
    parse_gff3 as pg_parse_gff3,
    parse_repeats,
    repeat_cds_coverage,
)


class TestCoveredBases:
    def test_no_overlap(self):
        assert covered_bases([(100, 200)], 300, 400) == 0

    def test_full_overlap(self):
        assert covered_bases([(1, 1000)], 100, 200) == 101

    def test_partial_overlap(self):
        assert covered_bases([(150, 250)], 100, 200) == 51

    def test_multiple_intervals(self):
        ivs = [(10, 50), (80, 120)]
        assert covered_bases(ivs, 1, 200) == 41 + 41


class TestRepeatCdsCoverage:
    def _make_children(self, cds_ranges):
        from detect_pseudogenes import GffRecord as R
        recs = [
            R('c', '.', 'CDS', s, e, '.', '+', '0', 'Parent=t1', '')
            for s, e in cds_ranges
        ]
        return {'t1': recs}

    def test_fully_covered(self):
        children = self._make_children([(100, 200)])
        cov = repeat_cds_coverage('t1', children, [(50, 300)])
        assert cov == pytest.approx(1.0)

    def test_not_covered(self):
        children = self._make_children([(100, 200)])
        cov = repeat_cds_coverage('t1', children, [(300, 400)])
        assert cov == pytest.approx(0.0)

    def test_partial_coverage(self):
        children = self._make_children([(100, 199)])  # 100 bp CDS
        cov = repeat_cds_coverage('t1', children, [(100, 149)])  # 50 bp overlap
        assert cov == pytest.approx(0.5)


class TestAllIntronsFrameshifted:
    def _make_children(self, exon_ranges):
        from detect_pseudogenes import GffRecord as R
        recs = [
            R('c', '.', 'exon', s, e, '.', '+', '.', 'Parent=t1', '')
            for s, e in exon_ranges
        ]
        return {'t1': recs}

    def test_single_exon_not_frameshifted(self):
        children = self._make_children([(1, 1000)])
        assert not all_introns_frameshifted('t1', children)

    def test_large_introns_not_frameshifted(self):
        # Intron = 499 bp (> FRAMESHIFT_INTRON_MAX=60)
        children = self._make_children([(1, 500), (1000, 2000)])
        assert not all_introns_frameshifted('t1', children)

    def test_all_small_introns_frameshifted(self):
        # Introns: 5 bp and 3 bp (both < 60)
        children = self._make_children([(1, 100), (106, 200), (204, 300)])
        assert all_introns_frameshifted('t1', children)

    def test_mixed_introns_not_all_frameshifted(self):
        # Introns: 5 bp and 200 bp
        children = self._make_children([(1, 100), (106, 200), (401, 500)])
        assert not all_introns_frameshifted('t1', children)


class TestClassifyGene:
    def _make_children_full(self, exon_ranges, cds_ranges):
        from detect_pseudogenes import GffRecord as R
        recs = []
        for s, e in exon_ranges:
            recs.append(R('c', '.', 'exon', s, e, '.', '+', '.', 'Parent=t1', ''))
        for s, e in cds_ranges:
            recs.append(R('c', '.', 'CDS', s, e, '.', '+', '0', 'Parent=t1', ''))
        return {'t1': recs}

    def test_single_exon_high_repeat_coverage(self):
        children = self._make_children_full([(100, 200)], [(100, 200)])
        result = classify_gene('g1', ['t1'], children, [(50, 300)], 0.80)
        assert result == 'processed_pseudogene'

    def test_single_exon_low_repeat_coverage(self):
        # Only 10% repeat coverage
        children = self._make_children_full([(100, 199)], [(100, 199)])  # 100 bp CDS
        result = classify_gene('g1', ['t1'], children, [(100, 109)], 0.80)  # 10 bp
        assert result is None

    def test_multi_exon_all_frameshifted_introns(self):
        # Two exons with 5 bp intron
        children = self._make_children_full([(1, 100), (106, 200)], [(1, 100), (106, 200)])
        result = classify_gene('g1', ['t1'], children, [], 0.80)
        assert result == 'pseudogene'

    def test_multi_exon_large_intron_not_pseudogene(self):
        # Two exons with 499 bp intron
        children = self._make_children_full([(1, 100), (600, 700)], [(1, 100), (600, 700)])
        result = classify_gene('g1', ['t1'], children, [], 0.80)
        assert result is None


# ===========================================================================
# ── detect_readthrough ──────────────────────────────────────────────────────
# ===========================================================================

from detect_readthrough import (
    GffRecord as RGffRecord,
    build_cds_index,
    is_readthrough,
    overlapping_genes,
    parse_gff3 as rt_parse_gff3,
)


def make_rt_gff3(tmp_path: Path, content: str) -> Path:
    p = tmp_path / 'input.gff3'
    p.write_text(textwrap.dedent(content))
    return p


class TestOverlappingGenes:
    def test_no_overlap(self):
        result = overlapping_genes(500, 600, [(100, 200, 'g1')], 'g0')
        assert result == set()

    def test_overlap_with_other_gene(self):
        result = overlapping_genes(100, 300, [(200, 400, 'g1')], 'g0')
        assert result == {'g1'}

    def test_overlap_with_own_gene_excluded(self):
        result = overlapping_genes(100, 300, [(200, 400, 'g1')], 'g1')
        assert result == set()

    def test_overlap_multiple_genes(self):
        cds = [(100, 200, 'g1'), (150, 250, 'g2')]
        result = overlapping_genes(50, 300, cds, 'g0')
        assert result == {'g1', 'g2'}


class TestIsReadthrough:
    def test_transcript_within_one_gene_not_readthrough(self, tmp_path):
        gff = make_rt_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t1\t1000\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t1\t1000\t.\t+\t0\tParent=t1
            chr1\t.\tgene\t2000\t3000\t.\t+\t.\tID=g2
            chr1\t.\tmRNA\t2000\t3000\t.\t+\t.\tID=t2;Parent=g2
            chr1\t.\texon\t2000\t3000\t.\t+\t.\tParent=t2
            chr1\t.\tCDS\t2000\t3000\t.\t+\t0\tParent=t2
        """)
        _, genes, transcripts, children, gene_order, tx_order = rt_parse_gff3(str(gff))
        cds_index = build_cds_index(gene_order, genes, tx_order, children)
        # t1 only overlaps g1's CDS region, not g2
        assert not is_readthrough('t1', 'g1', transcripts['t1'], children, cds_index)

    def test_transcript_spanning_two_genes_is_readthrough(self, tmp_path):
        gff = make_rt_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t4000\t.\t+\t.\tID=grt
            chr1\t.\tmRNA\t1\t4000\t.\t+\t.\tID=trt;Parent=grt
            chr1\t.\texon\t1\t4000\t.\t+\t.\tParent=trt
            chr1\t.\tCDS\t1\t4000\t.\t+\t0\tParent=trt
            chr1\t.\tgene\t500\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t500\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t500\t1000\t.\t+\t.\tParent=t1
            chr1\t.\tCDS\t500\t1000\t.\t+\t0\tParent=t1
            chr1\t.\tgene\t2000\t3000\t.\t+\t.\tID=g2
            chr1\t.\tmRNA\t2000\t3000\t.\t+\t.\tID=t2;Parent=g2
            chr1\t.\texon\t2000\t3000\t.\t+\t.\tParent=t2
            chr1\t.\tCDS\t2000\t3000\t.\t+\t0\tParent=t2
        """)
        _, genes, transcripts, children, gene_order, tx_order = rt_parse_gff3(str(gff))
        cds_index = build_cds_index(gene_order, genes, tx_order, children)
        # trt spans g1 and g2 CDSs
        assert is_readthrough('trt', 'grt', transcripts['trt'], children, cds_index)


# ===========================================================================
# ── select_canonical ────────────────────────────────────────────────────────
# ===========================================================================

from select_canonical import (
    GffRecord as CGffRecord,
    cds_length as sc_cds_length,
    parse_gff3 as sc_parse_gff3,
    select_canonical_tx,
    transcript_span,
)


def make_sc_gff3(tmp_path: Path, content: str) -> Path:
    p = tmp_path / 'input.gff3'
    p.write_text(textwrap.dedent(content))
    return p


class TestSelectCanonical:
    def test_longer_cds_wins(self, tmp_path):
        gff = make_sc_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t5000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t3000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t1\t900\t.\t+\t0\tParent=t1
            chr1\t.\tmRNA\t1\t5000\t.\t+\t.\tID=t2;Parent=g1
            chr1\t.\tCDS\t1\t1800\t.\t+\t0\tParent=t2
        """)
        _, _, transcripts, children, _, tx_order = sc_parse_gff3(str(gff))
        canonical = select_canonical_tx(tx_order['g1'], transcripts, children)
        assert canonical == 't2'

    def test_same_cds_longer_span_wins(self, tmp_path):
        # Both transcripts have same CDS length; t2 has a longer span
        gff = make_sc_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t5000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t100\t999\t.\t+\t0\tParent=t1
            chr1\t.\tmRNA\t1\t5000\t.\t+\t.\tID=t2;Parent=g1
            chr1\t.\tCDS\t100\t999\t.\t+\t0\tParent=t2
        """)
        _, _, transcripts, children, _, tx_order = sc_parse_gff3(str(gff))
        canonical = select_canonical_tx(tx_order['g1'], transcripts, children)
        assert canonical == 't2'

    def test_canonical_attribute_set_correctly(self, tmp_path):
        import subprocess
        gff_path = make_sc_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t5000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t3000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t1\t300\t.\t+\t0\tParent=t1
            chr1\t.\tmRNA\t1\t5000\t.\t+\t.\tID=t2;Parent=g1
            chr1\t.\tCDS\t1\t900\t.\t+\t0\tParent=t2
        """)
        out = tmp_path / 'out.gff3'
        script = os.path.join(BIN, 'select_canonical.py')
        result = subprocess.run(
            [sys.executable, script, '--gff3', str(gff_path), '--out', str(out)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        content = out.read_text()
        # t2 has longer CDS → canonical_transcript=1
        assert 'ID=t2' in content
        # find the t2 line and verify canonical_transcript=1
        for line in content.splitlines():
            if 'ID=t2' in line:
                assert 'canonical_transcript=1' in line
            if 'ID=t1' in line:
                assert 'canonical_transcript=0' in line

    def test_single_transcript_gene(self, tmp_path):
        gff = make_sc_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t1\t900\t.\t+\t0\tParent=t1
        """)
        _, _, transcripts, children, _, tx_order = sc_parse_gff3(str(gff))
        canonical = select_canonical_tx(tx_order['g1'], transcripts, children)
        assert canonical == 't1'

    def test_file_order_tiebreak(self, tmp_path):
        # Same CDS and same span — first in file (t1) should win
        gff = make_sc_gff3(tmp_path, """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\tCDS\t1\t300\t.\t+\t0\tParent=t1
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t2;Parent=g1
            chr1\t.\tCDS\t1\t300\t.\t+\t0\tParent=t2
        """)
        _, _, transcripts, children, _, tx_order = sc_parse_gff3(str(gff))
        canonical = select_canonical_tx(tx_order['g1'], transcripts, children)
        assert canonical == 't1'


# ===========================================================================
# ── flag_selenoproteins ─────────────────────────────────────────────────────
# ===========================================================================

from flag_selenoproteins import (
    parse_gff3 as sel_parse_gff3,
    parse_selenoprotein_names,
)


def make_sel_gff3(tmp_path: Path, content: str) -> Path:
    p = tmp_path / 'input.gff3'
    p.write_text(textwrap.dedent(content))
    return p


class TestFlagSelenoproteins:
    def _run_script(self, tmp_path, gff_content, fasta_content=None):
        import subprocess
        gff_path = make_sel_gff3(tmp_path, gff_content)
        if fasta_content is not None:
            fasta_path = tmp_path / 'seleno.fasta'
            fasta_path.write_text(textwrap.dedent(fasta_content))
        else:
            fasta_path = tmp_path / 'NO_FILE'
        out_path = tmp_path / 'out.gff3'
        script = os.path.join(BIN, 'flag_selenoproteins.py')
        result = subprocess.run(
            [sys.executable, script,
             '--gff3', str(gff_path),
             '--proteins', str(fasta_path),
             '--out', str(out_path)],
            capture_output=True, text=True,
        )
        return result, out_path

    def test_matching_gene_flagged(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=SELENOP;biotype=protein_coding
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1;biotype=protein_coding
        """
        fasta = """\
            >SELENOP selenoprotein P
            MKTIIALSYIFCLVFA
        """
        result, out_path = self._run_script(tmp_path, gff, fasta)
        assert result.returncode == 0, result.stderr
        content = out_path.read_text()
        assert 'biotype=selenoprotein' in content

    def test_non_matching_gene_unchanged(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=MYOGENE;biotype=protein_coding
            chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t1;Parent=g1;biotype=protein_coding
        """
        fasta = """\
            >SELENOP selenoprotein P
            MKTIIALSYIFCLVFA
        """
        result, out_path = self._run_script(tmp_path, gff, fasta)
        assert result.returncode == 0, result.stderr
        content = out_path.read_text()
        assert 'biotype=selenoprotein' not in content
        assert 'biotype=protein_coding' in content

    def test_no_file_passthrough(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=g1;Name=SELENOP;biotype=protein_coding
        """
        # fasta_content=None → uses NO_FILE sentinel
        result, out_path = self._run_script(tmp_path, gff, fasta_content=None)
        assert result.returncode == 0, result.stderr
        # Output should be identical to input
        original = (tmp_path / 'input.gff3').read_text()
        assert out_path.read_text() == original

    def test_parse_selenoprotein_names(self, tmp_path):
        fasta = tmp_path / 'seleno.fasta'
        fasta.write_text(textwrap.dedent("""\
            >SELENOP selenoprotein P description here
            MKTIIALSYIFCLVFA
            >GPX4 glutathione peroxidase 4
            ACGT
        """))
        names = parse_selenoprotein_names(str(fasta))
        assert 'SELENOP' in names
        assert 'GPX4' in names
        assert len(names) == 2


# ===========================================================================
# ── Integration: filter_geneset end-to-end via subprocess ───────────────────
# ===========================================================================

class TestFilterGenesetIntegration:
    def _run(self, tmp_path, gff_content, min_orf_aa=100, min_intron_size=10):
        import subprocess
        gff_path = tmp_path / 'input.gff3'
        gff_path.write_text(textwrap.dedent(gff_content))
        out_path = tmp_path / 'out.filtered.gff3'
        script = os.path.join(BIN, 'filter_geneset.py')
        result = subprocess.run(
            [sys.executable, script,
             '--gff3', str(gff_path),
             '--out', str(out_path),
             '--min-orf-aa', str(min_orf_aa),
             '--min-intron-size', str(min_intron_size)],
            capture_output=True, text=True,
        )
        return result, out_path

    def test_short_orf_transcript_removed(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t_short;Parent=g1
            chr1\t.\texon\t1\t2000\t.\t+\t.\tParent=t_short
            chr1\t.\tCDS\t1\t60\t.\t+\t0\tParent=t_short
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t_long;Parent=g1
            chr1\t.\texon\t1\t2000\t.\t+\t.\tParent=t_long
            chr1\t.\tCDS\t1\t900\t.\t+\t0\tParent=t_long
        """
        result, out_path = self._run(tmp_path, gff, min_orf_aa=100)
        assert result.returncode == 0, result.stderr
        content = out_path.read_text()
        assert 't_short' not in content
        assert 't_long' in content

    def test_tiny_intron_transcript_removed(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t2000\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t_tiny;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t_tiny
            chr1\t.\texon\t504\t2000\t.\t+\t.\tParent=t_tiny
            chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t_tiny
            chr1\t.\tCDS\t504\t2000\t.\t+\t0\tParent=t_tiny
            chr1\t.\tmRNA\t1\t2000\t.\t+\t.\tID=t_ok;Parent=g1
            chr1\t.\texon\t1\t500\t.\t+\t.\tParent=t_ok
            chr1\t.\texon\t600\t2000\t.\t+\t.\tParent=t_ok
            chr1\t.\tCDS\t1\t500\t.\t+\t0\tParent=t_ok
            chr1\t.\tCDS\t600\t2000\t.\t+\t0\tParent=t_ok
        """
        result, out_path = self._run(tmp_path, gff, min_intron_size=10)
        assert result.returncode == 0, result.stderr
        content = out_path.read_text()
        assert 't_tiny' not in content
        assert 't_ok' in content

    def test_gene_removed_when_all_transcripts_filtered(self, tmp_path):
        gff = """\
            ##gff-version 3
            chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g_gone
            chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t_gone;Parent=g_gone
            chr1\t.\texon\t1\t100\t.\t+\t.\tParent=t_gone
            chr1\t.\tCDS\t1\t60\t.\t+\t0\tParent=t_gone
        """
        result, out_path = self._run(tmp_path, gff, min_orf_aa=100)
        assert result.returncode == 0, result.stderr
        content = out_path.read_text()
        assert 'g_gone' not in content
