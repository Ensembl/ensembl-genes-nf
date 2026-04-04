"""
Tests for project_transcripts.py
"""

import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from project_transcripts import (
    ChainBlock,
    Chain,
    SourceTranscript,
    parse_chain,
    parse_source_gff3,
    _project_interval,
    _compute_coverage,
    project,
    write_projected_gff3,
)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

CHAIN_SIMPLE = textwrap.dedent("""\
    chain 1000 chr1 100000 + 0 100000 chr1 100000 + 0 100000 1
    10000

""")

# Chain with gap: two blocks separated by dt=500 (target gap) and dq=0 (no query gap)
CHAIN_WITH_GAP = textwrap.dedent("""\
    chain 800 chr1 100000 + 1000 12000 chr2 100000 + 1000 12000 2
    5000\t500\t0
    5000

""")

# Minus-strand chain: query maps to - strand
CHAIN_MINUS = textwrap.dedent("""\
    chain 900 chr1 100000 + 0 10000 chr2 100000 - 90000 100000 3
    10000

""")

GFF3_BASIC = textwrap.dedent("""\
    ##gff-version 3
    chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;biotype=protein_coding
    chr1\tEnsembl\ttranscript\t1000\t5000\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding
    chr1\tEnsembl\texon\t1000\t2000\t.\t+\t.\tID=tx1_exon1;Parent=tx1
    chr1\tEnsembl\texon\t3000\t5000\t.\t+\t.\tID=tx1_exon2;Parent=tx1
    chr1\tEnsembl\tgene\t8000\t10000\t.\t-\t.\tID=gene2;biotype=lncRNA
    chr1\tEnsembl\ttranscript\t8000\t10000\t.\t-\t.\tID=tx2;Parent=gene2;biotype=lncRNA
    chr1\tEnsembl\texon\t8000\t10000\t.\t-\t.\tID=tx2_exon1;Parent=tx2
""")


def _write_tmp(content: str, suffix='.chain') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _simple_chain(q_strand='+', score=1000, t_start=0, t_end=10000,
                  q_start=0, q_end=10000, q_size=100000) -> Chain:
    c = Chain(
        score=score,
        t_name='chr1', t_size=100000, t_start=t_start, t_end=t_end,
        q_name='chr2', q_size=q_size, q_strand=q_strand,
        q_start=q_start, q_end=q_end,
        chain_id='test',
    )
    block_size = t_end - t_start
    if q_strand == '+':
        q_block_start = q_start
    else:
        q_block_start = q_size - q_end
    c.blocks.append(ChainBlock(
        t_start=t_start, t_end=t_end,
        q_start=q_block_start, q_end=q_block_start + block_size,
        q_strand=q_strand,
    ))
    return c


# ---------------------------------------------------------------------------
# parse_chain tests
# ---------------------------------------------------------------------------

class TestParseChain:
    def test_parses_single_chain(self):
        path = _write_tmp(CHAIN_SIMPLE)
        chains = parse_chain(path)
        os.unlink(path)
        assert len(chains) == 1

    def test_chain_header_fields(self):
        path = _write_tmp(CHAIN_SIMPLE)
        chains = parse_chain(path)
        os.unlink(path)
        c = chains[0]
        assert c.score == 1000
        assert c.t_name == 'chr1'
        assert c.q_name == 'chr1'
        assert c.t_start == 0
        assert c.t_end == 100000

    def test_single_block_parsed(self):
        path = _write_tmp(CHAIN_SIMPLE)
        chains = parse_chain(path)
        os.unlink(path)
        assert len(chains[0].blocks) == 1
        b = chains[0].blocks[0]
        assert b.t_start == 0
        assert b.t_end == 10000

    def test_chain_with_gap_has_two_blocks(self):
        path = _write_tmp(CHAIN_WITH_GAP)
        chains = parse_chain(path)
        os.unlink(path)
        assert len(chains[0].blocks) == 2

    def test_gap_advances_t_pos(self):
        path = _write_tmp(CHAIN_WITH_GAP)
        chains = parse_chain(path)
        os.unlink(path)
        b1, b2 = chains[0].blocks
        # b1: t_start=1000, t_end=6000; gap dt=500; b2: t_start=6500
        assert b1.t_start == 1000
        assert b1.t_end == 6000
        assert b2.t_start == 6500
        assert b2.t_end == 11500

    def test_minus_strand_chain(self):
        path = _write_tmp(CHAIN_MINUS)
        chains = parse_chain(path)
        os.unlink(path)
        assert chains[0].q_strand == '-'

    def test_minus_strand_block_coords(self):
        path = _write_tmp(CHAIN_MINUS)
        chains = parse_chain(path)
        os.unlink(path)
        b = chains[0].blocks[0]
        # Chain: q_start=90000, q_end=100000, q_size=100000, block_size=10000
        # parse_chain sets q_pos = q_size - q_end = 0
        # q_block_start = q_size - q_pos - block_size = 100000 - 0 - 10000 = 90000
        # Block covers [90000, 100000) in forward query coords
        assert b.q_start == 90000
        assert b.q_end == 100000

    def test_empty_file(self):
        path = _write_tmp('')
        chains = parse_chain(path)
        os.unlink(path)
        assert chains == []


# ---------------------------------------------------------------------------
# parse_source_gff3 tests
# ---------------------------------------------------------------------------

class TestParseSourceGff3:
    def test_parses_transcripts(self):
        path = _write_tmp(GFF3_BASIC, suffix='.gff3')
        txs = parse_source_gff3(path)
        os.unlink(path)
        assert len(txs) == 2

    def test_exons_attached(self):
        path = _write_tmp(GFF3_BASIC, suffix='.gff3')
        txs = parse_source_gff3(path)
        os.unlink(path)
        tx = next(t for t in txs if t.tx_id == 'tx1')
        assert len(tx.exons) == 2
        assert tx.exons[0] == (1000, 2000)
        assert tx.exons[1] == (3000, 5000)

    def test_gene_id_from_parent(self):
        path = _write_tmp(GFF3_BASIC, suffix='.gff3')
        txs = parse_source_gff3(path)
        os.unlink(path)
        tx = next(t for t in txs if t.tx_id == 'tx1')
        assert tx.gene_id == 'gene1'

    def test_strand_parsed(self):
        path = _write_tmp(GFF3_BASIC, suffix='.gff3')
        txs = parse_source_gff3(path)
        os.unlink(path)
        tx2 = next(t for t in txs if t.tx_id == 'tx2')
        assert tx2.strand == '-'

    def test_empty_gff3(self):
        path = _write_tmp('##gff-version 3\n', suffix='.gff3')
        txs = parse_source_gff3(path)
        os.unlink(path)
        assert txs == []


# ---------------------------------------------------------------------------
# _project_interval tests
# ---------------------------------------------------------------------------

class TestProjectInterval:
    def _chain(self, q_strand='+'):
        return _simple_chain(q_strand=q_strand, t_start=0, t_end=10000,
                              q_start=0, q_end=10000)

    def test_full_overlap(self):
        c = self._chain()
        result = _project_interval(1, 10000, c)
        assert result is not None
        q_s, q_e, strand = result
        assert q_s == 1
        assert q_e == 10000
        assert strand == '+'

    def test_partial_overlap(self):
        c = self._chain()
        # t_start=1000, t_end=2000 (1-based)
        result = _project_interval(1001, 2000, c)
        assert result is not None
        q_s, q_e, strand = result
        assert q_s == 1001
        assert q_e == 2000

    def test_no_overlap_returns_none(self):
        c = self._chain()
        # interval is outside chain range
        result = _project_interval(50001, 60000, c)
        assert result is None

    def test_minus_strand_strand_returned(self):
        c = self._chain(q_strand='-')
        result = _project_interval(1, 10000, c)
        assert result is not None
        assert result[2] == '-'

    def test_1based_conversion(self):
        # interval [1, 1000] maps to [1, 1000] in a 1:1 chain starting at 0
        c = self._chain()
        q_s, q_e, _ = _project_interval(1, 1000, c)
        assert q_s == 1
        assert q_e == 1000


# ---------------------------------------------------------------------------
# _compute_coverage tests
# ---------------------------------------------------------------------------

class TestComputeCoverage:
    def test_all_projected(self):
        exons = [(1000, 2000), (3000, 5000)]
        proj  = [(1000, 2000, '+'), (3000, 5000, '+')]
        cov = _compute_coverage(exons, proj)
        assert cov == pytest.approx(100.0)

    def test_one_exon_missing(self):
        exons = [(1000, 2000), (3000, 5000)]   # 1001 + 2001 = 3002 bases
        proj  = [(1000, 2000, '+'), None]        # only 1001 bases mapped
        cov = _compute_coverage(exons, proj)
        assert cov == pytest.approx(1001 / 3002 * 100.0)

    def test_all_missing_returns_zero(self):
        exons = [(1000, 2000)]
        proj  = [None]
        cov = _compute_coverage(exons, proj)
        assert cov == 0.0

    def test_empty_exons(self):
        cov = _compute_coverage([], [])
        assert cov == 0.0


# ---------------------------------------------------------------------------
# project (integration) tests
# ---------------------------------------------------------------------------

class TestProject:
    def _source_tx(self, seqname='chr1', start=1001, end=5000,
                   exons=None, tx_id='tx1', gene_id='gene1'):
        tx = SourceTranscript(
            seqname=seqname, start=start, end=end,
            strand='+', tx_id=tx_id, gene_id=gene_id,
            biotype='protein_coding',
        )
        tx.exons = exons or [(1001, 2000), (3001, 5000)]
        return tx

    def test_projects_transcript_above_threshold(self):
        chain = _simple_chain(t_start=0, t_end=10000, q_start=0, q_end=10000)
        txs = [self._source_tx()]
        results = project(txs, [chain], min_coverage=50.0)
        assert len(results) == 1

    def test_filters_below_coverage_threshold(self):
        # Only one exon overlaps chain (covers < 50%)
        # chain covers [0, 2000), tx exons are [1001,2000] + [3001,5000]
        chain = _simple_chain(t_start=0, t_end=2000, q_start=0, q_end=2000)
        txs = [self._source_tx()]
        results = project(txs, [chain], min_coverage=50.0)
        # exon1: 1000 bp, exon2: 2000 bp → total 3000; only 1000 maps → ~33%
        assert len(results) == 0

    def test_no_chain_for_seqname(self):
        chain = _simple_chain(t_start=0, t_end=10000)
        chain.t_name = 'chrX'
        tx = self._source_tx(seqname='chr1')
        results = project([tx], [chain], min_coverage=0.0)
        assert results == []

    def test_picks_best_chain_by_score(self):
        c1 = _simple_chain(score=500,  t_start=0, t_end=10000, q_start=0, q_end=10000)
        c2 = _simple_chain(score=2000, t_start=0, t_end=10000, q_start=20000, q_end=30000)
        c1.q_name = 'chrA'
        c2.q_name = 'chrB'
        txs = [self._source_tx()]
        results = project(txs, [c1, c2], min_coverage=0.0)
        assert len(results) == 1
        assert results[0]['q_name'] == 'chrB'

    def test_result_fields_present(self):
        chain = _simple_chain(t_start=0, t_end=10000, q_start=0, q_end=10000)
        txs = [self._source_tx()]
        result = project(txs, [chain], min_coverage=0.0)[0]
        for key in ('src_tx_id', 'src_gene_id', 'q_name', 'q_start',
                    'q_end', 'q_strand', 'coverage', 'proj_exons'):
            assert key in result

    def test_coverage_100_when_all_exons_map(self):
        chain = _simple_chain(t_start=0, t_end=10000, q_start=0, q_end=10000)
        txs = [self._source_tx(exons=[(1001, 2000), (3001, 5000)])]
        result = project(txs, [chain], min_coverage=0.0)[0]
        assert result['coverage'] == pytest.approx(100.0)


# ---------------------------------------------------------------------------
# write_projected_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteProjectedGff3:
    def _projections(self):
        chain = _simple_chain(t_start=0, t_end=10000, q_start=0, q_end=10000)
        txs = [
            SourceTranscript(
                seqname='chr1', start=1001, end=5000,
                strand='+', tx_id='tx1', gene_id='gene1',
                biotype='protein_coding',
                exons=[(1001, 2000), (3001, 5000)],
            ),
            SourceTranscript(
                seqname='chr1', start=6001, end=8000,
                strand='+', tx_id='tx2', gene_id='gene1',
                biotype='protein_coding',
                exons=[(6001, 8000)],
            ),
        ]
        return project(txs, [chain], min_coverage=0.0)

    def test_writes_header(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            assert fh.readline().strip() == '##gff-version 3'
        os.unlink(out.name)

    def test_biotype_projected_transcript(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert 'biotype=projected_transcript' in content
        os.unlink(out.name)

    def test_gene_feature_written(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            lines = fh.readlines()
        assert any('\tgene\t' in l for l in lines)
        os.unlink(out.name)

    def test_exon_features_written(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            content = fh.read()
        assert '\texon\t' in content
        os.unlink(out.name)

    def test_returns_tx_count(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        n = write_projected_gff3(projs, out.name)
        assert n == len(projs)
        os.unlink(out.name)

    def test_src_gene_grouped(self):
        # Both txs share gene1 → should produce one gene feature
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        assert len(gene_lines) == 1
        os.unlink(out.name)

    def test_gene_id_prefix(self):
        projs = self._projections()
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        write_projected_gff3(projs, out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        for line in gene_lines:
            gene_id = line.split('ID=')[1].split(';')[0]
            assert gene_id.startswith('proj_gene_')
        os.unlink(out.name)
