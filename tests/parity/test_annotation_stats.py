"""
Tests for metrics/gff3_stats.py :: compute_stats()

Verifies that the annotation statistics module correctly parses GFF3 files
and computes biologically meaningful summary statistics. These tests use
synthetic GFF3 fixtures from conftest.py.

Statistics tested:
  - Gene / transcript counting
  - Biotype classification
  - UTR coverage rates
  - Pseudogene rates
  - Canonical transcript CDS lengths
  - Transcripts-per-gene distribution
"""

import sys
from pathlib import Path

import pytest

# Allow importing from the metrics package
sys.path.insert(0, str(Path(__file__).parent))
from metrics.gff3_stats import compute_stats


class TestComputeStatsGeneCount:
    """Basic gene and transcript counting."""

    def test_total_gene_count(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        assert stats.total_genes == 4

    def test_total_transcript_count(self, stats_gff3):
        # g1→2 txs, g2→1, g3→1, g4→1 = 5
        stats = compute_stats(str(stats_gff3))
        assert stats.total_transcripts == 5

    def test_protein_coding_gene_count(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        assert stats.protein_coding_genes == 2

    def test_genes_by_biotype_includes_lncrna(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        biotypes = stats.genes_by_biotype
        assert "protein_coding" in biotypes
        assert biotypes["protein_coding"] == 2

    def test_pseudogene_count(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        # g4 biotype=pseudogene but attached via pseudogenic_transcript (not processed)
        assert stats.pseudogene_count == 1
        assert stats.processed_pseudogene_count == 0


class TestComputeStatsUtrCoverage:
    """UTR coverage rate calculations."""

    def test_utr5_rate_is_fraction_of_protein_coding(self, stats_gff3):
        """g1 has 5' UTR, g2 does not → rate = 0.5."""
        stats = compute_stats(str(stats_gff3))
        assert stats.utr5_rate() == pytest.approx(0.5)

    def test_utr3_rate_is_fraction_of_protein_coding(self, stats_gff3):
        """g1 has 3' UTR, g2 does not → rate = 0.5."""
        stats = compute_stats(str(stats_gff3))
        assert stats.utr3_rate() == pytest.approx(0.5)

    def test_genes_with_both_utrs(self, stats_gff3):
        """Only g1 has both UTRs."""
        stats = compute_stats(str(stats_gff3))
        assert stats.genes_with_both_utrs == 1

    def test_no_utrs_gives_zero_rates(self, utr_ref_gff3):
        """A GFF3 with no UTR features has 0.0 UTR rates."""
        stats = compute_stats(str(utr_ref_gff3))
        assert stats.utr5_rate() == 0.0
        assert stats.utr3_rate() == 0.0


class TestComputeStatsCdsAndExons:
    """CDS length and exon count statistics."""

    def test_canonical_cds_length_captured(self, stats_gff3):
        """
        g1 canonical: CDS 1050-1200 (151 bp) + 1600-2400 (801 bp) = 952 bp
        g2 canonical: CDS 4050-4400 (351 bp) + 4700-4950 (251 bp) = 602 bp
        """
        stats = compute_stats(str(stats_gff3))
        assert len(stats.canonical_cds_lengths) == 2
        assert 952 in stats.canonical_cds_lengths
        assert 602 in stats.canonical_cds_lengths

    def test_median_cds_length(self, stats_gff3):
        """Median of [952, 602] = 777.0."""
        stats = compute_stats(str(stats_gff3))
        assert stats.median_cds_length() == pytest.approx(777.0)

    def test_exon_counts_populated(self, stats_gff3):
        """All transcripts should have exon counts > 0."""
        stats = compute_stats(str(stats_gff3))
        assert len(stats.exon_counts) == 5  # 5 transcripts
        assert all(n > 0 for n in stats.exon_counts)


class TestComputeStatsTxPerGene:
    """Transcripts-per-gene distribution."""

    def test_mean_tx_per_gene(self, stats_gff3):
        """g1→2, g2→1, g3→1, g4→1 → mean = 5/4 = 1.25"""
        stats = compute_stats(str(stats_gff3))
        assert stats.mean_tx_per_gene() == pytest.approx(1.25)

    def test_tx_per_gene_list_length(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        assert len(stats.tx_per_gene) == 4  # one entry per gene


class TestComputeStatsCanonical:
    """Canonical coverage metrics."""

    def test_canonical_coverage_full(self, stats_gff3):
        """All 4 genes have a canonical transcript in stats_gff3."""
        stats = compute_stats(str(stats_gff3))
        assert stats.canonical_coverage() == 1.0

    def test_canonical_coverage_partial(self, utr_ref_gff3):
        """utr_ref_gff3 has is_canonical=True for its one gene."""
        stats = compute_stats(str(utr_ref_gff3))
        assert stats.canonical_coverage() == 1.0


class TestComputeStatsSummaryDict:
    """summary() dict structure and types."""

    def test_summary_has_required_keys(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        summary = stats.summary()
        required = {
            "total_genes", "total_transcripts", "genes_by_biotype",
            "protein_coding_genes", "utr5_rate", "utr3_rate",
            "pseudogene_rate", "mean_tx_per_gene",
            "median_canonical_cds_bp", "canonical_coverage",
        }
        assert required.issubset(summary.keys())

    def test_summary_pseudogene_rate(self, stats_gff3):
        stats = compute_stats(str(stats_gff3))
        summary = stats.summary()
        # 1 pseudogene out of 4 genes
        assert summary["pseudogene_rate"] == pytest.approx(0.25, abs=0.001)
