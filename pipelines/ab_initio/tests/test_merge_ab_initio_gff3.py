"""
Tests for merge_ab_initio_gff3.py
"""

import os
import re
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from merge_ab_initio_gff3 import (
    _update_ids,
    count_genes,
    merge,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GFF3_CHUNK1 = textwrap.dedent("""\
    ##gff-version 3
    chr1\taugustus\tgene\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000001;biotype=ab_initio
    chr1\taugustus\ttranscript\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000001_tx_0001;Parent=ab_initio_gene_00000001
    chr1\taugustus\texon\t1000\t2000\t.\t+\t.\tID=ab_initio_gene_00000001_tx_0001_exon_1;Parent=ab_initio_gene_00000001_tx_0001
""")

GFF3_CHUNK2 = textwrap.dedent("""\
    ##gff-version 3
    chr2\taugustus\tgene\t100\t500\t.\t-\t.\tID=ab_initio_gene_00000001;biotype=ab_initio
    chr2\taugustus\ttranscript\t100\t500\t.\t-\t.\tID=ab_initio_gene_00000001_tx_0001;Parent=ab_initio_gene_00000001
    chr2\taugustus\texon\t100\t500\t.\t-\t.\tID=ab_initio_gene_00000001_tx_0001_exon_1;Parent=ab_initio_gene_00000001_tx_0001
    chr2\taugustus\tgene\t600\t900\t.\t+\t.\tID=ab_initio_gene_00000002;biotype=ab_initio
    chr2\taugustus\ttranscript\t600\t900\t.\t+\t.\tID=ab_initio_gene_00000002_tx_0001;Parent=ab_initio_gene_00000002
    chr2\taugustus\texon\t600\t900\t.\t+\t.\tID=ab_initio_gene_00000002_tx_0001_exon_1;Parent=ab_initio_gene_00000002_tx_0001
""")


def _write_tmp(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.gff3', delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


# ---------------------------------------------------------------------------
# _update_ids tests
# ---------------------------------------------------------------------------

class TestUpdateIds:
    def test_gene_id_incremented(self):
        lines = ['chr1\taugustus\tgene\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000001;biotype=ab_initio']
        result = _update_ids(lines, gene_offset=5)
        assert 'ab_initio_gene_00000006' in result[0]

    def test_parent_updated_consistently(self):
        lines = [
            'chr1\taugustus\tgene\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000001',
            'chr1\taugustus\ttranscript\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000001_tx_0001;Parent=ab_initio_gene_00000001',
        ]
        result = _update_ids(lines, gene_offset=2)
        assert 'ID=ab_initio_gene_00000003' in result[0]
        assert 'Parent=ab_initio_gene_00000003' in result[1]

    def test_zero_offset_unchanged(self):
        lines = ['chr1\taugustus\tgene\t1000\t5000\t.\t+\t.\tID=ab_initio_gene_00000003']
        result = _update_ids(lines, gene_offset=0)
        assert result[0] == lines[0]

    def test_comment_lines_preserved(self):
        lines = ['# comment', 'chr1\taugustus\tgene\t1\t100\t.\t+\t.\tID=ab_initio_gene_00000001']
        result = _update_ids(lines, gene_offset=0)
        assert result[0] == '# comment'


# ---------------------------------------------------------------------------
# count_genes tests
# ---------------------------------------------------------------------------

class TestCountGenes:
    def test_counts_gene_lines(self):
        lines = [
            'chr1\taugustus\tgene\t1\t100\t.\t+\t.\tID=ab_initio_gene_00000001',
            'chr1\taugustus\ttranscript\t1\t100\t.\t+\t.\tID=tx1',
            'chr2\taugustus\tgene\t200\t300\t.\t-\t.\tID=ab_initio_gene_00000002',
        ]
        assert count_genes(lines) == 2

    def test_ignores_comments(self):
        lines = ['# gene line not real', 'chr1\taugustus\tgene\t1\t100\t.\t+\t.\tID=g1']
        assert count_genes(lines) == 1

    def test_empty(self):
        assert count_genes([]) == 0


# ---------------------------------------------------------------------------
# merge (integration) tests
# ---------------------------------------------------------------------------

class TestMerge:
    def test_writes_header(self):
        p1 = _write_tmp(GFF3_CHUNK1)
        p2 = _write_tmp(GFF3_CHUNK2)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        merge([p1, p2], out.name)
        with open(out.name) as fh:
            assert fh.readline().strip() == '##gff-version 3'
        os.unlink(p1); os.unlink(p2); os.unlink(out.name)

    def test_returns_total_gene_count(self):
        p1 = _write_tmp(GFF3_CHUNK1)
        p2 = _write_tmp(GFF3_CHUNK2)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        n = merge([p1, p2], out.name)
        # chunk1 has 1 gene, chunk2 has 2 genes
        assert n == 3
        os.unlink(p1); os.unlink(p2); os.unlink(out.name)

    def test_ids_sequential(self):
        p1 = _write_tmp(GFF3_CHUNK1)
        p2 = _write_tmp(GFF3_CHUNK2)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        merge([p1, p2], out.name)
        with open(out.name) as fh:
            gene_lines = [l for l in fh if '\tgene\t' in l]
        ids = [l.split('ID=')[1].split(';')[0] for l in gene_lines]
        # Should be 00000001, 00000002, 00000003 (chunk2's genes get offset by 1)
        assert ids[0] == 'ab_initio_gene_00000001'
        assert ids[1] == 'ab_initio_gene_00000002'
        assert ids[2] == 'ab_initio_gene_00000003'
        os.unlink(p1); os.unlink(p2); os.unlink(out.name)

    def test_parent_refs_consistent(self):
        p1 = _write_tmp(GFF3_CHUNK1)
        p2 = _write_tmp(GFF3_CHUNK2)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        merge([p1, p2], out.name)
        with open(out.name) as fh:
            lines = [l.rstrip() for l in fh if '\t' in l]

        # For each non-gene line, verify its Parent exists as an ID somewhere
        ids = set()
        for line in lines:
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            attrs = cols[8]
            m = re.search(r'ID=([^;]+)', attrs)
            if m:
                ids.add(m.group(1))

        for line in lines:
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            m = re.search(r'Parent=([^;]+)', cols[8])
            if m:
                assert m.group(1) in ids, f"Parent {m.group(1)} not found in IDs"

        os.unlink(p1); os.unlink(p2); os.unlink(out.name)

    def test_single_chunk(self):
        p1 = _write_tmp(GFF3_CHUNK1)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        n = merge([p1], out.name)
        assert n == 1
        os.unlink(p1); os.unlink(out.name)
