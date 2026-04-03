"""Unit tests for bin/chunk_fasta.py"""

import sys
import os
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))
from chunk_fasta import parse_fasta, write_chunk


class TestParseFasta(unittest.TestCase):

    def _write_fasta(self, entries):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        for seqid, seq in entries:
            f.write(f">{seqid}\n{seq}\n")
        f.close()
        return f.name

    def test_single_sequence(self):
        tmp = self._write_fasta([('chr1', 'ATCGATCG')])
        result = list(parse_fasta(tmp))
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][0], 'chr1')
        self.assertEqual(result[0][1], 'ATCGATCG')
        os.unlink(tmp)

    def test_multiple_sequences(self):
        tmp = self._write_fasta([('chr1', 'AAAA'), ('chr2', 'CCCC'), ('chr3', 'GGGG')])
        result = list(parse_fasta(tmp))
        self.assertEqual(len(result), 3)
        self.assertEqual(result[1][0], 'chr2')
        os.unlink(tmp)

    def test_header_with_description_uses_first_word(self):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        f.write(">chr1 some description here\nATCG\n")
        f.close()
        result = list(parse_fasta(f.name))
        self.assertEqual(result[0][0], 'chr1')
        os.unlink(f.name)

    def test_multiline_sequence_joined(self):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        f.write(">chr1\nAAAA\nCCCC\nGGGG\n")
        f.close()
        result = list(parse_fasta(f.name))
        self.assertEqual(result[0][1], 'AAAACCCCGGGG')
        os.unlink(f.name)


class TestWriteChunk(unittest.TestCase):

    def test_chunk_file_written(self):
        with tempfile.TemporaryDirectory() as d:
            out = write_chunk('chr1', 0, 100, 'A' * 100, Path(d))
            self.assertTrue(os.path.exists(out))

    def test_chunk_header_contains_coordinates(self):
        with tempfile.TemporaryDirectory() as d:
            out = write_chunk('chr1', 100, 200, 'T' * 100, Path(d))
            with open(out) as f:
                header = f.readline().strip()
            self.assertEqual(header, '>chr1_100_200')

    def test_chunk_sequence_correct_length(self):
        with tempfile.TemporaryDirectory() as d:
            seq = 'ACGT' * 25  # 100 bp
            out = write_chunk('chr1', 0, 100, seq, Path(d))
            with open(out) as f:
                lines = [l.strip() for l in f if not l.startswith('>')]
            actual = ''.join(lines)
            self.assertEqual(actual, seq)

    def test_partial_chunk_at_end(self):
        with tempfile.TemporaryDirectory() as d:
            seq = 'N' * 150
            out = write_chunk('chr1', 100, 150, seq, Path(d))
            with open(out) as f:
                lines = [l.strip() for l in f if not l.startswith('>')]
            actual = ''.join(lines)
            self.assertEqual(len(actual), 50)


class TestMergeRepeats(unittest.TestCase):
    """Test merge_repeats.py parsing helpers."""

    def setUp(self):
        sys.path.insert(0, str(Path(__file__).parents[1] / 'bin'))

    def test_source_from_filename(self):
        from merge_repeats import source_from_filename
        self.assertEqual(source_from_filename('genome.rpt.bed'), 'RED')
        self.assertEqual(source_from_filename('genome.trf.bed'), 'TRF')
        self.assertEqual(source_from_filename('genome.dust.bed'), 'DUST')
        self.assertEqual(source_from_filename('genome.fa.out'), 'RepeatMasker')

    def test_parse_bed_line_valid(self):
        from merge_repeats import parse_bed_line
        result = parse_bed_line('chr1\t100\t200\tLINE\t0\t+\n', 'RepeatMasker')
        self.assertIsNotNone(result)
        self.assertEqual(result[0], 'chr1')
        self.assertEqual(result[1], 100)
        self.assertEqual(result[2], 200)

    def test_parse_bed_line_too_short(self):
        from merge_repeats import parse_bed_line
        result = parse_bed_line('chr1\t100\n', 'RepeatMasker')
        self.assertIsNone(result)

    def test_bed_to_gff3_format(self):
        from merge_repeats import bed_to_gff3
        bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed.write('chr1\t100\t200\tLINE\t0\t.\tRepeatMasker\n')
        bed.write('chr1\t500\t600\tSINE\t0\t.\tRepeatMasker\n')
        bed.close()
        gff = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        gff.close()
        n = bed_to_gff3(bed.name, gff.name, 'test')
        self.assertEqual(n, 2)
        with open(gff.name) as f:
            lines = [l for l in f if not l.startswith('#')]
        cols = lines[0].split('\t')
        self.assertEqual(cols[2], 'repeat_region')
        self.assertEqual(int(cols[3]), 101)  # GFF3 1-based
        os.unlink(bed.name)
        os.unlink(gff.name)


if __name__ == '__main__':
    unittest.main()
