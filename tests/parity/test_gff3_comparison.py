"""
Tests for metrics/compare.py :: compare_gff3()

Verifies that the GFF3 comparison module correctly:
  1. Matches genes by CDS overlap (not gene name)
  2. Computes CDS boundary accuracy
  3. Detects UTR additions/removals
  4. Scores pseudogene recall and precision
  5. Counts novel and missed genes
  6. Handles identical annotations correctly (all metrics at 1.0)

Uses synthetic GFF3 fixtures from conftest.py, each encoding a specific
comparison scenario with known expected results.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from metrics.compare import compare_gff3, ComparisonResult


class TestIdenticalAnnotations:
    """When test == ref (same CDS, different IDs), all metrics should be perfect."""

    def test_all_ref_genes_matched(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        assert result.total_ref_genes == 2
        assert result.matched_genes == 2

    def test_no_unmatched_or_novel_genes(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        assert result.unmatched_ref_genes == 0
        assert result.novel_test_genes == 0

    def test_cds_boundary_accuracy_is_perfect(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        assert result.cds_boundary_accuracy == pytest.approx(1.0)
        assert result.exact_cds_match == 2

    def test_biotype_fully_concordant(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        assert result.biotype_concordant == 2
        assert result.biotype_discordant == 0

    def test_utr_deltas_are_zero(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        assert result.utr5_added == [0, 0]
        assert result.utr3_added == [0, 0]

    def test_gene_recall_and_precision_are_one(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        summary = result.summary()
        assert summary["gene_recall"] == pytest.approx(1.0)
        assert summary["gene_precision"] == pytest.approx(1.0)


class TestUtrAddition:
    """New pipeline adds UTR — CDS should still match, UTR delta should be positive."""

    def test_gene_matched_by_cds_overlap(self, utr_ref_gff3, utr_test_gff3):
        result = compare_gff3(str(utr_ref_gff3), str(utr_test_gff3))
        assert result.matched_genes == 1
        assert result.unmatched_ref_genes == 0

    def test_cds_boundary_accuracy_still_perfect(self, utr_ref_gff3, utr_test_gff3):
        """UTR addition does not change CDS → CDS boundaries should still match."""
        result = compare_gff3(str(utr_ref_gff3), str(utr_test_gff3))
        assert result.cds_boundary_accuracy == pytest.approx(1.0)

    def test_utr5_positive_delta(self, utr_ref_gff3, utr_test_gff3):
        """Test pipeline added 150 bp of 5' UTR that ref lacks."""
        result = compare_gff3(str(utr_ref_gff3), str(utr_test_gff3))
        assert len(result.utr5_added) == 1
        assert result.utr5_added[0] == 150   # 1000-1049 = 150 bp

    def test_utr3_positive_delta(self, utr_ref_gff3, utr_test_gff3):
        """Test pipeline added 150 bp of 3' UTR that ref lacks."""
        result = compare_gff3(str(utr_ref_gff3), str(utr_test_gff3))
        assert len(result.utr3_added) == 1
        assert result.utr3_added[0] == 150   # 1951-2100 = 150 bp

    def test_summary_mean_utr_delta_positive(self, utr_ref_gff3, utr_test_gff3):
        result = compare_gff3(str(utr_ref_gff3), str(utr_test_gff3))
        summary = result.summary()
        assert summary["mean_utr5_delta_bp"] == 150.0
        assert summary["mean_utr3_delta_bp"] == 150.0


class TestCdsBoundaryShift:
    """New pipeline has a shifted CDS boundary — should detect as non-exact CDS match."""

    def test_gene_matched_despite_cds_shift(self, shifted_ref_gff3, shifted_test_gff3):
        """CDS overlap still holds even with 3 bp shift."""
        result = compare_gff3(str(shifted_ref_gff3), str(shifted_test_gff3))
        assert result.matched_genes == 1

    def test_cds_boundary_accuracy_is_zero(self, shifted_ref_gff3, shifted_test_gff3):
        """Shifted CDS → no exact match."""
        result = compare_gff3(str(shifted_ref_gff3), str(shifted_test_gff3))
        assert result.exact_cds_match == 0
        assert result.partial_cds_match == 1
        assert result.cds_boundary_accuracy == pytest.approx(0.0)

    def test_summary_shows_low_cds_accuracy(self, shifted_ref_gff3, shifted_test_gff3):
        result = compare_gff3(str(shifted_ref_gff3), str(shifted_test_gff3))
        summary = result.summary()
        assert summary["cds_boundary_accuracy"] == pytest.approx(0.0)


class TestNovelGene:
    """New pipeline calls a gene with no reference counterpart."""

    def test_matched_ref_gene_found(self, novel_gene_ref_gff3, novel_gene_test_gff3):
        result = compare_gff3(str(novel_gene_ref_gff3), str(novel_gene_test_gff3))
        assert result.matched_genes == 1

    def test_novel_gene_counted(self, novel_gene_ref_gff3, novel_gene_test_gff3):
        """The chr2 gene in test has no match in ref."""
        result = compare_gff3(str(novel_gene_ref_gff3), str(novel_gene_test_gff3))
        assert result.novel_test_genes == 1

    def test_precision_below_one(self, novel_gene_ref_gff3, novel_gene_test_gff3):
        result = compare_gff3(str(novel_gene_ref_gff3), str(novel_gene_test_gff3))
        summary = result.summary()
        # 1 matched out of 2 test genes
        assert summary["gene_precision"] == pytest.approx(0.5)

    def test_recall_still_one(self, novel_gene_ref_gff3, novel_gene_test_gff3):
        result = compare_gff3(str(novel_gene_ref_gff3), str(novel_gene_test_gff3))
        summary = result.summary()
        assert summary["gene_recall"] == pytest.approx(1.0)


class TestMissedGene:
    """New pipeline misses a gene present in reference (false negative)."""

    def test_matched_one_of_two_ref_genes(self, missed_gene_ref_gff3, missed_gene_test_gff3):
        result = compare_gff3(str(missed_gene_ref_gff3), str(missed_gene_test_gff3))
        assert result.matched_genes == 1

    def test_unmatched_ref_gene_counted(self, missed_gene_ref_gff3, missed_gene_test_gff3):
        result = compare_gff3(str(missed_gene_ref_gff3), str(missed_gene_test_gff3))
        assert result.unmatched_ref_genes == 1

    def test_no_novel_test_genes(self, missed_gene_ref_gff3, missed_gene_test_gff3):
        result = compare_gff3(str(missed_gene_ref_gff3), str(missed_gene_test_gff3))
        assert result.novel_test_genes == 0

    def test_recall_below_one(self, missed_gene_ref_gff3, missed_gene_test_gff3):
        result = compare_gff3(str(missed_gene_ref_gff3), str(missed_gene_test_gff3))
        summary = result.summary()
        assert summary["gene_recall"] == pytest.approx(0.5)


class TestPseudogeneRecall:
    """Pseudogene detection — recall and precision scoring."""

    def test_ref_pseudogene_count(self, pseudo_ref_gff3, pseudo_missed_test_gff3):
        result = compare_gff3(str(pseudo_ref_gff3), str(pseudo_missed_test_gff3))
        assert result.ref_pseudogenes == 1

    def test_pseudogene_not_recalled(self, pseudo_ref_gff3, pseudo_missed_test_gff3):
        """Test pipeline called protein_coding instead of pseudogene → recall=0."""
        result = compare_gff3(str(pseudo_ref_gff3), str(pseudo_missed_test_gff3))
        assert result.pseudogene_recall == pytest.approx(0.0)

    def test_biotype_discordant_for_pseudogene_locus(self, pseudo_ref_gff3, pseudo_missed_test_gff3):
        """The pseudo locus is matched by CDS but biotypes differ."""
        result = compare_gff3(str(pseudo_ref_gff3), str(pseudo_missed_test_gff3))
        # gene1 matches (pc→pc concordant), pseudo1 matches (pseudo→pc discordant)
        assert result.biotype_discordant >= 1

    def test_test_pseudogene_count_is_zero(self, pseudo_ref_gff3, pseudo_missed_test_gff3):
        result = compare_gff3(str(pseudo_ref_gff3), str(pseudo_missed_test_gff3))
        assert result.test_pseudogenes == 0


class TestComparisonResultSummary:
    """summary() dict structure, types, and edge cases."""

    def test_summary_keys_complete(self, identical_ref_gff3, identical_test_gff3):
        result = compare_gff3(str(identical_ref_gff3), str(identical_test_gff3))
        summary = result.summary()
        required = {
            "ref_genes", "test_genes", "matched_genes", "gene_recall",
            "gene_precision", "unmatched_ref_genes", "novel_test_genes",
            "cds_boundary_accuracy", "biotype_concordance",
            "mean_utr5_delta_bp", "mean_utr3_delta_bp",
            "pseudogene_recall", "pseudogene_precision",
        }
        assert required.issubset(summary.keys())

    def test_summary_utr_delta_none_when_no_matches(self, tmp_path):
        """If no genes matched, UTR delta lists are empty → mean should be None."""
        # ref has one gene, test has a gene on different chromosome → no match
        ref = tmp_path / "r.gff3"
        ref.write_text(
            "##gff-version 3\n"
            "chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=g1;biotype=protein_coding\n"
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
            "chr1\t.\tCDS\t1050\t1950\t.\t+\t0\tID=c1;Parent=t1\n"
        )
        test = tmp_path / "t.gff3"
        test.write_text(
            "##gff-version 3\n"
            "chrX\t.\tgene\t1000\t2000\t.\t+\t.\tID=g2;biotype=protein_coding\n"
            "chrX\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t2;Parent=g2\n"
            "chrX\t.\tCDS\t1050\t1950\t.\t+\t0\tID=c2;Parent=t2\n"
        )
        result = compare_gff3(str(ref), str(test))
        summary = result.summary()
        assert summary["mean_utr5_delta_bp"] is None
        assert summary["mean_utr3_delta_bp"] is None
        assert result.matched_genes == 0
