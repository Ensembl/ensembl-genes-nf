"""
Unit tests for bin/collapse_long_reads.py

Tests the pure-Python logic (BAM parsing helpers, collapse, GFF3 output)
without requiring pysam or real BAM files — pysam is mocked where needed.
"""

import sys
import os
import tempfile
import unittest
from pathlib import Path

# Add bin/ to path so we can import the script directly
sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))

from collapse_long_reads import (
    Exon,
    Transcript,
    cluster_by_locus,
    collapse_cluster,
    write_gff3,
    write_stats,
    _cigar_to_exons,
)


# ---------------------------------------------------------------------------
# Exon
# ---------------------------------------------------------------------------

class TestExon(unittest.TestCase):

    def test_length(self):
        e = Exon(100, 200)
        self.assertEqual(e.length(), 100)

    def test_overlap_full(self):
        e1 = Exon(100, 200)
        e2 = Exon(100, 200)
        self.assertEqual(e1.overlap(e2), 100)

    def test_overlap_partial(self):
        e1 = Exon(100, 200)
        e2 = Exon(150, 250)
        self.assertEqual(e1.overlap(e2), 50)

    def test_overlap_none(self):
        e1 = Exon(100, 200)
        e2 = Exon(300, 400)
        self.assertEqual(e1.overlap(e2), 0)


# ---------------------------------------------------------------------------
# Transcript
# ---------------------------------------------------------------------------

class TestTranscript(unittest.TestCase):

    def _make(self, chrom, strand, exon_coords):
        exons = [Exon(s, e) for s, e in exon_coords]
        return Transcript('t1', chrom, strand, exons)

    def test_start_end(self):
        t = self._make('chr1', '+', [(100, 200), (300, 400)])
        self.assertEqual(t.start, 100)
        self.assertEqual(t.end,   400)

    def test_exon_length(self):
        t = self._make('chr1', '+', [(100, 200), (300, 400)])
        self.assertEqual(t.exon_length(), 200)

    def test_reciprocal_overlap_identical(self):
        t1 = self._make('chr1', '+', [(100, 200), (300, 400)])
        t2 = self._make('chr1', '+', [(100, 200), (300, 400)])
        self.assertAlmostEqual(t1.reciprocal_overlap(t2), 1.0)

    def test_reciprocal_overlap_none(self):
        t1 = self._make('chr1', '+', [(100, 200)])
        t2 = self._make('chr1', '+', [(300, 400)])
        self.assertAlmostEqual(t1.reciprocal_overlap(t2), 0.0)

    def test_reciprocal_overlap_different_strand(self):
        t1 = self._make('chr1', '+', [(100, 200)])
        t2 = self._make('chr1', '-', [(100, 200)])
        self.assertAlmostEqual(t1.reciprocal_overlap(t2), 0.0)

    def test_reciprocal_overlap_different_chrom(self):
        t1 = self._make('chr1', '+', [(100, 200)])
        t2 = self._make('chr2', '+', [(100, 200)])
        self.assertAlmostEqual(t1.reciprocal_overlap(t2), 0.0)

    def test_reciprocal_overlap_partial(self):
        # t1 exon: 100-300 (length 200), t2 exon: 200-400 (length 200)
        # shared: 200-300 = 100
        # reciprocal = 100 / min(200, 200) = 0.5
        t1 = self._make('chr1', '+', [(100, 300)])
        t2 = self._make('chr1', '+', [(200, 400)])
        self.assertAlmostEqual(t1.reciprocal_overlap(t2), 0.5)


# ---------------------------------------------------------------------------
# _cigar_to_exons
# ---------------------------------------------------------------------------

class TestCigarToExons(unittest.TestCase):
    # CIGAR ops: 0=M, 3=N(intron)
    M, N = 0, 3

    def test_simple_match(self):
        exons = _cigar_to_exons(100, [(self.M, 200)], max_intron=200000)
        self.assertEqual(len(exons), 1)
        self.assertEqual(exons[0].start, 100)
        self.assertEqual(exons[0].end,   300)

    def test_spliced(self):
        # 100M 1000N 100M starting at ref pos 0
        exons = _cigar_to_exons(0, [(self.M, 100), (self.N, 1000), (self.M, 100)], max_intron=200000)
        self.assertEqual(len(exons), 2)
        self.assertEqual(exons[0].start, 0)
        self.assertEqual(exons[0].end,   100)
        self.assertEqual(exons[1].start, 1100)
        self.assertEqual(exons[1].end,   1200)

    def test_intron_too_large_returns_empty(self):
        exons = _cigar_to_exons(0, [(self.M, 100), (self.N, 500000), (self.M, 100)], max_intron=200000)
        self.assertEqual(exons, [])


# ---------------------------------------------------------------------------
# cluster_by_locus
# ---------------------------------------------------------------------------

class TestClusterByLocus(unittest.TestCase):

    def _tx(self, chrom, strand, start, end, tid='t'):
        return Transcript(tid, chrom, strand, [Exon(start, end)])

    def test_two_overlapping(self):
        t1 = self._tx('chr1', '+', 100, 500, 't1')
        t2 = self._tx('chr1', '+', 300, 700, 't2')
        clusters = cluster_by_locus([t1, t2])
        self.assertEqual(len(clusters), 1)
        self.assertEqual(len(clusters[0]), 2)

    def test_two_non_overlapping(self):
        t1 = self._tx('chr1', '+', 100, 200, 't1')
        t2 = self._tx('chr1', '+', 500, 600, 't2')
        clusters = cluster_by_locus([t1, t2])
        self.assertEqual(len(clusters), 2)

    def test_different_strands_separate(self):
        t1 = self._tx('chr1', '+', 100, 500, 't1')
        t2 = self._tx('chr1', '-', 100, 500, 't2')
        clusters = cluster_by_locus([t1, t2])
        self.assertEqual(len(clusters), 2)

    def test_different_chroms_separate(self):
        t1 = self._tx('chr1', '+', 100, 500, 't1')
        t2 = self._tx('chr2', '+', 100, 500, 't2')
        clusters = cluster_by_locus([t1, t2])
        self.assertEqual(len(clusters), 2)


# ---------------------------------------------------------------------------
# collapse_cluster
# ---------------------------------------------------------------------------

class TestCollapseCluster(unittest.TestCase):

    def _tx(self, tid, start, end):
        return Transcript(tid, 'chr1', '+', [Exon(start, end)])

    def test_identical_transcripts_collapse_to_one(self):
        txs = [self._tx(f't{i}', 100, 300) for i in range(5)]
        result = collapse_cluster(txs, min_overlap=0.8)
        self.assertEqual(len(result), 1)

    def test_non_overlapping_all_kept(self):
        t1 = self._tx('t1', 100, 200)
        t2 = self._tx('t2', 500, 600)
        result = collapse_cluster([t1, t2], min_overlap=0.8)
        self.assertEqual(len(result), 2)

    def test_partial_overlap_below_threshold_kept(self):
        # 50% overlap, threshold 0.8 -> both kept
        t1 = self._tx('t1', 100, 300)  # length 200
        t2 = self._tx('t2', 200, 400)  # length 200, 100bp overlap -> 0.5
        result = collapse_cluster([t1, t2], min_overlap=0.8)
        self.assertEqual(len(result), 2)

    def test_partial_overlap_above_threshold_collapsed(self):
        # 90% overlap, threshold 0.8 -> merged
        t1 = Transcript('t1', 'chr1', '+', [Exon(100, 200)])  # length 100
        t2 = Transcript('t2', 'chr1', '+', [Exon(100, 192)])  # length 92, 92/100=0.92
        result = collapse_cluster([t1, t2], min_overlap=0.8)
        self.assertEqual(len(result), 1)

    def test_empty_cluster(self):
        result = collapse_cluster([], min_overlap=0.8)
        self.assertEqual(result, [])


# ---------------------------------------------------------------------------
# write_gff3
# ---------------------------------------------------------------------------

class TestWriteGff3(unittest.TestCase):

    def _tx(self, tid, chrom, strand, exons):
        return Transcript(tid, chrom, strand, [Exon(s, e) for s, e in exons])

    def test_output_is_valid_gff3(self):
        txs = [
            self._tx('t1', 'chr1', '+', [(100, 200), (300, 400)]),
            self._tx('t2', 'chr1', '-', [(500, 600)]),
        ]
        with tempfile.NamedTemporaryFile(mode='r', suffix='.gff3', delete=False) as f:
            tmp = f.name

        write_gff3(txs, tmp, 'test_sample')

        with open(tmp) as f:
            lines = [l for l in f if not l.startswith('#') and l.strip()]

        feature_types = [l.split('\t')[2] for l in lines]
        self.assertIn('gene',       feature_types)
        self.assertIn('transcript', feature_types)
        self.assertIn('exon',       feature_types)
        os.unlink(tmp)

    def test_gene_ids_are_unique(self):
        txs = [
            self._tx(f't{i}', 'chr1', '+', [(i * 1000, i * 1000 + 100)])
            for i in range(1, 4)
        ]
        with tempfile.NamedTemporaryFile(mode='r', suffix='.gff3', delete=False) as f:
            tmp = f.name

        write_gff3(txs, tmp, 'sample')
        with open(tmp) as f:
            gene_lines = [l for l in f if '\tgene\t' in l]

        gene_ids = [l.split('ID=')[1].split(';')[0] for l in gene_lines]
        self.assertEqual(len(gene_ids), len(set(gene_ids)))
        os.unlink(tmp)


# ---------------------------------------------------------------------------
# write_stats
# ---------------------------------------------------------------------------

class TestWriteStats(unittest.TestCase):

    def test_stats_format(self):
        txs = [Transcript(f't{i}', 'chr1', '+', [Exon(i*100, i*100+50)]) for i in range(3)]
        with tempfile.NamedTemporaryFile(mode='r', suffix='.tsv', delete=False) as f:
            tmp = f.name

        write_stats(txs, tmp, 'mysample')
        with open(tmp) as f:
            lines = f.readlines()

        self.assertEqual(lines[0].strip(), 'sample\tgenes\ttranscripts')
        cols = lines[1].strip().split('\t')
        self.assertEqual(cols[0], 'mysample')
        self.assertTrue(cols[1].isdigit())
        self.assertTrue(cols[2].isdigit())
        os.unlink(tmp)


if __name__ == '__main__':
    unittest.main()
