"""
Tests for filter_exonerate.py
"""

import gzip
import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from filter_exonerate import (
    ExonerateHit,
    parse_exonerate_gff,
    load_query_lengths,
    write_gff3,
    _compute_stats,
    _open,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gff(lines: str, suffix='.gff') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(lines)
    tmp.close()
    return tmp.name


def _make_fasta(seqs: dict, suffix='.fa') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    for sid, seq in seqs.items():
        tmp.write(f'>{sid}\n{seq}\n')
    tmp.close()
    return tmp.name


GFF_CDNA = textwrap.dedent("""\
    ##gff-version 2
    # exonerate output
    chr1\texonerate:cdna2genome\tgene\t1000\t5000\t500\t+\t.\tQuery NM_001 ; score 500
    chr1\texonerate:cdna2genome\texon\t1000\t2000\t.\t+\t.\tsequence NM_001
    chr1\texonerate:cdna2genome\texon\t3000\t5000\t.\t+\t.\tsequence NM_001
    chr1\texonerate:cdna2genome\tgene\t6000\t8000\t300\t-\t.\tQuery NM_002 ; score 300
    chr1\texonerate:cdna2genome\texon\t7000\t8000\t.\t-\t.\tsequence NM_002
    chr1\texonerate:cdna2genome\texon\t6000\t6500\t.\t-\t.\tsequence NM_002
""")

GFF_PROTEIN = textwrap.dedent("""\
    ##gff-version 2
    chr1\texonerate:protein2genome\tgene\t1200\t4800\t300\t+\t.\tQuery P12345 ; score 300
    chr1\texonerate:protein2genome\texon\t1200\t2200\t.\t+\t.\tsequence P12345
    chr1\texonerate:protein2genome\texon\t3000\t4800\t.\t+\t.\tsequence P12345
""")

GFF_MINUS_STRAND = textwrap.dedent("""\
    chr1\texonerate:cdna2genome\tgene\t5000\t1000\t400\t-\t.\tQuery NM_003 ; score 400
    chr1\texonerate:cdna2genome\texon\t5000\t4000\t.\t-\t.\tsequence NM_003
""")


# ---------------------------------------------------------------------------
# parse_exonerate_gff tests
# ---------------------------------------------------------------------------

class TestParseExonerateGff:
    def test_parses_cdna_genes(self):
        path = _make_gff(GFF_CDNA)
        hits = parse_exonerate_gff(path, 'cdna')
        assert len(hits) == 2
        os.unlink(path)

    def test_cdna_hit_coordinates(self):
        path = _make_gff(GFF_CDNA)
        hits = parse_exonerate_gff(path, 'cdna')
        h = hits[0]
        assert h.seqname == 'chr1'
        assert h.start == 1000
        assert h.end   == 5000
        assert h.strand == '+'
        assert h.score == 500.0
        assert h.query_id == 'NM_001'
        os.unlink(path)

    def test_cdna_hit_exons(self):
        path = _make_gff(GFF_CDNA)
        hits = parse_exonerate_gff(path, 'cdna')
        h = hits[0]
        assert len(h.exons) == 2
        assert h.exons[0] == (1000, 2000, '+')
        assert h.exons[1] == (3000, 5000, '+')
        os.unlink(path)

    def test_protein_biotype(self):
        path = _make_gff(GFF_PROTEIN)
        hits = parse_exonerate_gff(path, 'protein')
        assert hits[0].biotype == 'protein_alignment'
        os.unlink(path)

    def test_cdna_biotype(self):
        path = _make_gff(GFF_CDNA)
        hits = parse_exonerate_gff(path, 'cdna')
        assert hits[0].biotype == 'cdna_alignment'
        os.unlink(path)

    def test_minus_strand_coord_normalisation(self):
        path = _make_gff(GFF_MINUS_STRAND)
        hits = parse_exonerate_gff(path, 'cdna')
        h = hits[0]
        assert h.start <= h.end, "start must be <= end after normalisation"
        assert h.start == 1000
        assert h.end   == 5000
        os.unlink(path)

    def test_empty_file_returns_empty(self):
        path = _make_gff('##gff-version 2\n')
        hits = parse_exonerate_gff(path, 'cdna')
        assert hits == []
        os.unlink(path)

    def test_comment_lines_ignored(self):
        gff = '# this is a comment\n' + GFF_CDNA
        path = _make_gff(gff)
        hits = parse_exonerate_gff(path, 'cdna')
        assert len(hits) == 2
        os.unlink(path)

    def test_gz_support(self):
        path = _make_gff('', suffix='.gff.gz')
        with gzip.open(path, 'wt') as fh:
            fh.write(GFF_CDNA)
        hits = parse_exonerate_gff(path, 'cdna')
        assert len(hits) == 2
        os.unlink(path)


# ---------------------------------------------------------------------------
# load_query_lengths tests
# ---------------------------------------------------------------------------

class TestLoadQueryLengths:
    def test_basic_lengths(self):
        path = _make_fasta({'SEQ1': 'ATCG' * 10, 'SEQ2': 'ATCG' * 5})
        lengths = load_query_lengths(path)
        assert lengths['SEQ1'] == 40
        assert lengths['SEQ2'] == 20
        os.unlink(path)

    def test_none_returns_empty(self):
        assert load_query_lengths(None) == {}

    def test_multiline_fasta(self):
        path = _make_fasta({'SEQ1': 'ATCG'})
        # append a multiline manually
        with open(path, 'a') as fh:
            fh.write('>SEQ2\nATCG\nATCG\n')
        lengths = load_query_lengths(path)
        assert lengths['SEQ2'] == 8
        os.unlink(path)


# ---------------------------------------------------------------------------
# _compute_stats tests
# ---------------------------------------------------------------------------

class TestComputeStats:
    def _hit(self, exons, score=500.0, biotype='cdna_alignment', query_id='NM_001'):
        return ExonerateHit(
            seqname='chr1', start=1000, end=5000, strand='+',
            score=score, coverage=0.0, pid=0.0,
            query_id=query_id, biotype=biotype, exons=exons,
        )

    def test_coverage_with_known_query_length(self):
        hit = self._hit([(1000, 2000, '+'), (3000, 5000, '+')])
        # exon1 = 2000-1000+1 = 1001 bp, exon2 = 5000-3000+1 = 2001 bp → total 3002 bp
        query_lengths = {'NM_001': 3002}
        _compute_stats(hit, query_lengths)
        assert abs(hit.coverage - 100.0) < 1.0

    def test_coverage_without_known_length(self):
        hit = self._hit([(1000, 2000, '+'), (3000, 5000, '+')])
        _compute_stats(hit, {})
        assert 0.0 < hit.coverage <= 100.0

    def test_pid_above_zero(self):
        hit = self._hit([(1000, 2000, '+')], score=500)
        _compute_stats(hit, {})
        assert hit.pid > 0.0

    def test_coverage_capped_at_100(self):
        hit = self._hit([(1000, 2000, '+')], query_id='NM_SHORT')
        _compute_stats(hit, {'NM_SHORT': 1})  # query shorter than exon
        assert hit.coverage <= 100.0

    def test_protein_biotype_uses_higher_pts(self):
        hit_cdna    = self._hit([(1000, 2000, '+')], score=500, biotype='cdna_alignment')
        hit_protein = self._hit([(1000, 2000, '+')], score=500, biotype='protein_alignment')
        _compute_stats(hit_cdna, {})
        _compute_stats(hit_protein, {})
        # protein uses 10 pts/match vs 5 for cdna → lower pid for same score
        assert hit_protein.pid <= hit_cdna.pid


# ---------------------------------------------------------------------------
# write_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def _make_hits(self):
        h = ExonerateHit(
            seqname='chr1', start=1000, end=5000, strand='+',
            score=500.0, coverage=95.0, pid=97.0,
            query_id='NM_001', biotype='cdna_alignment',
            exons=[(1000, 2000, '+'), (3000, 5000, '+')],
        )
        return [h]

    def test_writes_gff_version_header(self):
        tmp = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        tmp.close()
        write_gff3(self._make_hits(), tmp.name, 'CHUNK01')
        with open(tmp.name) as fh:
            content = fh.read()
        assert '##gff-version 3' in content
        os.unlink(tmp.name)

    def test_writes_gene_transcript_exon(self):
        tmp = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        tmp.close()
        write_gff3(self._make_hits(), tmp.name, 'CHUNK01')
        with open(tmp.name) as fh:
            lines = [l for l in fh if '\t' in l]
        features = [l.split('\t')[2] for l in lines]
        assert 'gene' in features
        assert 'transcript' in features
        assert 'exon' in features
        os.unlink(tmp.name)

    def test_returns_count(self):
        tmp = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        tmp.close()
        n = write_gff3(self._make_hits(), tmp.name, 'CHUNK01')
        assert n == 1
        os.unlink(tmp.name)

    def test_empty_hits_writes_only_header(self):
        tmp = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        tmp.close()
        n = write_gff3([], tmp.name, 'CHUNK01')
        assert n == 0
        with open(tmp.name) as fh:
            content = fh.read().strip()
        assert content == '##gff-version 3'
        os.unlink(tmp.name)
