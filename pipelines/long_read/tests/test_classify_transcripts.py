"""
Unit tests for bin/classify_transcripts.py
"""

import sys
import os
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))

from classify_transcripts import (
    parse_blast,
    classify_gff3,
    _get_attr,
    _set_attr,
)


SAMPLE_GFF3 = """\
##gff-version 3
chr1\tlong_read\tgene\t101\t400\t.\t+\t.\tID=sample_gene_000001;biotype=isoseq
chr1\tlong_read\ttranscript\t101\t400\t.\t+\t.\tID=sample_gene_000001.t1;Parent=sample_gene_000001;biotype=isoseq
chr1\tlong_read\texon\t101\t200\t.\t+\t.\tParent=sample_gene_000001.t1
chr1\tlong_read\texon\t301\t400\t.\t+\t.\tParent=sample_gene_000001.t1
chr1\tlong_read\tgene\t1001\t1200\t.\t+\t.\tID=sample_gene_000002;biotype=isoseq
chr1\tlong_read\ttranscript\t1001\t1200\t.\t+\t.\tID=sample_gene_000002.t1;Parent=sample_gene_000002;biotype=isoseq
chr1\tlong_read\texon\t1001\t1200\t.\t+\t.\tParent=sample_gene_000002.t1
"""


class TestGetSetAttr(unittest.TestCase):

    def test_get_existing(self):
        self.assertEqual(_get_attr('ID=foo;biotype=isoseq', 'ID'), 'foo')
        self.assertEqual(_get_attr('ID=foo;biotype=isoseq', 'biotype'), 'isoseq')

    def test_get_missing(self):
        self.assertEqual(_get_attr('ID=foo', 'biotype'), '')

    def test_set_existing(self):
        result = _set_attr('ID=foo;biotype=isoseq', 'biotype', 'isoseq_supported')
        self.assertIn('biotype=isoseq_supported', result)
        self.assertNotIn('biotype=isoseq;', result)

    def test_set_new(self):
        result = _set_attr('ID=foo', 'biotype', 'isoseq')
        self.assertIn('biotype=isoseq', result)


class TestParseBlast(unittest.TestCase):

    def _write_blast(self, rows):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        for row in rows:
            f.write('\t'.join(str(c) for c in row) + '\n')
        f.close()
        return f.name

    def test_supported_ids_below_evalue(self):
        rows = [
            ['sample_gene_000001.t1', 'sp|P12345|GENE_HUMAN', 95.0, 100, 5, 0, 1, 100, 1, 100, '1e-10', 200],
        ]
        tmp = self._write_blast(rows)
        supported = parse_blast(tmp, evalue_cutoff=1e-5)
        self.assertIn('sample_gene_000001.t1', supported)
        os.unlink(tmp)

    def test_above_evalue_excluded(self):
        rows = [
            ['sample_gene_000001.t1', 'sp|P12345|GENE_HUMAN', 30.0, 50, 30, 5, 1, 50, 1, 50, '0.1', 20],
        ]
        tmp = self._write_blast(rows)
        supported = parse_blast(tmp, evalue_cutoff=1e-5)
        self.assertNotIn('sample_gene_000001.t1', supported)
        os.unlink(tmp)

    def test_empty_blast_file(self):
        tmp = self._write_blast([])
        supported = parse_blast(tmp, evalue_cutoff=1e-5)
        self.assertEqual(len(supported), 0)
        os.unlink(tmp)

    def test_comment_lines_skipped(self):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
        f.write('# this is a comment\n')
        f.write('tx1\tsp|X\t99\t100\t1\t0\t1\t100\t1\t100\t1e-50\t300\n')
        f.close()
        supported = parse_blast(f.name, evalue_cutoff=1e-5)
        self.assertIn('tx1', supported)
        os.unlink(f.name)


class TestClassifyGff3(unittest.TestCase):

    def setUp(self):
        self.gff3_in = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
        self.gff3_in.write(SAMPLE_GFF3)
        self.gff3_in.close()
        self.gff3_out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        self.gff3_out.close()

    def tearDown(self):
        os.unlink(self.gff3_in.name)
        os.unlink(self.gff3_out.name)

    def _read_out(self):
        with open(self.gff3_out.name) as f:
            return f.read()

    def test_supported_transcript_gets_biotype(self):
        supported = {'sample_gene_000001.t1'}
        classify_gff3(self.gff3_in.name, self.gff3_out.name, supported)
        content = self._read_out()
        # The supported transcript line should have isoseq_supported
        tx_lines = [l for l in content.splitlines() if '\ttranscript\t' in l]
        biotypes = [_get_attr(l.split('\t')[8], 'biotype') for l in tx_lines]
        self.assertIn('isoseq_supported', biotypes)

    def test_unsupported_transcript_keeps_isoseq(self):
        supported = set()  # nothing supported
        classify_gff3(self.gff3_in.name, self.gff3_out.name, supported)
        content = self._read_out()
        tx_lines = [l for l in content.splitlines() if '\ttranscript\t' in l]
        biotypes = [_get_attr(l.split('\t')[8], 'biotype') for l in tx_lines]
        self.assertTrue(all(b == 'isoseq' for b in biotypes))

    def test_gene_biotype_updated_when_child_supported(self):
        supported = {'sample_gene_000001.t1'}
        counts = classify_gff3(self.gff3_in.name, self.gff3_out.name, supported)
        self.assertEqual(counts['genes_supported'], 1)
        content = self._read_out()
        gene_lines = [l for l in content.splitlines() if '\tgene\t' in l]
        # gene_000001 should be isoseq_supported, gene_000002 should still be isoseq
        g1 = next(l for l in gene_lines if 'gene_000001' in l)
        g2 = next(l for l in gene_lines if 'gene_000002' in l)
        self.assertIn('isoseq_supported', _get_attr(g1.split('\t')[8], 'biotype'))
        self.assertEqual(_get_attr(g2.split('\t')[8], 'biotype'), 'isoseq')

    def test_counts_correct(self):
        supported = {'sample_gene_000001.t1'}
        counts = classify_gff3(self.gff3_in.name, self.gff3_out.name, supported)
        self.assertEqual(counts['isoseq_supported'], 1)
        self.assertEqual(counts['isoseq'], 1)
        self.assertEqual(counts['genes_supported'], 1)


if __name__ == '__main__':
    unittest.main()
