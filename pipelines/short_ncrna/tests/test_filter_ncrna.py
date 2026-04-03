"""Tests for filter_ncrna.py"""

from __future__ import annotations

import sys, os, textwrap
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from filter_ncrna import (
    rfam_biotype, CmsearchHit, parse_tblout, BlastHit,
    parse_blast_tabular, write_gff3,
)


# ---------------------------------------------------------------------------
# rfam_biotype
# ---------------------------------------------------------------------------

class TestRfamBiotype:
    def test_mirna(self):
        assert rfam_biotype('mir-21') == 'miRNA'

    def test_snorna(self):
        assert rfam_biotype('SNORD14') == 'snoRNA'

    def test_scarna_before_snorna(self):
        assert rfam_biotype('scaRNA9') == 'scaRNA'

    def test_snrna_u1(self):
        assert rfam_biotype('U1') == 'snRNA'

    def test_rrna_lsu(self):
        assert rfam_biotype('LSU_rRNA_eukarya') == 'rRNA'

    def test_rrna_5s(self):
        assert rfam_biotype('5S_rRNA') == 'rRNA'

    def test_trna(self):
        assert rfam_biotype('tRNA-Ala') == 'tRNA'

    def test_ribozyme(self):
        assert rfam_biotype('Hammerhead_1') == 'ribozyme'

    def test_vault(self):
        assert rfam_biotype('Vault') == 'vault_RNA'

    def test_misc_rna_fallback(self):
        assert rfam_biotype('UNKNOWN_RF99999') == 'misc_RNA'


# ---------------------------------------------------------------------------
# CmsearchHit.from_line
# ---------------------------------------------------------------------------

# Columns: target_name  tacc  query_name qacc mdl mdl_from mdl_to
#          seq_from seq_to strand trunc pass gc bias score evalue inc desc...
TBLOUT_GOOD = (
    'chr1\t-\tmir-21\tRF00XXX\tcm\t1\t88'
    '\t1000\t1088\t+\t-\t1\t0.46\t0.0\t87.3\t1.2e-16\t!\tmiRNA mir-21\n'
)
TBLOUT_MINUS = (
    'chr2\t-\tSNORD14\tRF00YYY\tcm\t1\t65'
    '\t5000\t4936\t-\t-\t1\t0.40\t0.0\t72.0\t3.0e-10\t!\tsnoRNA\n'
)
TBLOUT_BAD = '# this is a comment line\n'
TBLOUT_SHORT = 'chr1\t-\tmir-21\n'


class TestCmsearchHitFromLine:
    def test_good_line_parsed(self):
        hit = CmsearchHit.from_line(TBLOUT_GOOD)
        assert hit is not None
        assert hit.seqname == 'chr1'
        assert hit.cm_name == 'mir-21'
        assert hit.seq_from == 1000
        assert hit.seq_to == 1088
        assert hit.strand == '+'
        assert hit.evalue == pytest.approx(1.2e-16)
        assert hit.biotype == 'miRNA'

    def test_minus_strand_coords_normalised(self):
        hit = CmsearchHit.from_line(TBLOUT_MINUS)
        assert hit is not None
        assert hit.seq_from < hit.seq_to  # normalised
        assert hit.strand == '-'
        assert hit.biotype == 'snoRNA'

    def test_comment_returns_none(self):
        assert CmsearchHit.from_line(TBLOUT_BAD) is None

    def test_short_line_returns_none(self):
        assert CmsearchHit.from_line(TBLOUT_SHORT) is None


# ---------------------------------------------------------------------------
# parse_tblout (with tmp file)
# ---------------------------------------------------------------------------

TBLOUT_FILE = textwrap.dedent("""\
    # cmsearch :: search CM(s) against a sequence database
    chr1\t-\tmir-21\tRF001\tcm\t1\t88\t1000\t1088\t+\t-\t1\t0.46\t0.0\t87.3\t1.2e-16\t!\tmiRNA
    chr1\t-\tSNORD14\tRF002\tcm\t1\t65\t2000\t2065\t+\t-\t1\t0.40\t0.0\t72.0\t3.0e-10\t!\tsnoRNA
    chr1\t-\tLSU_rRNA\tRF003\tcm\t1\t200\t3000\t3200\t+\t-\t1\t0.45\t0.0\t1750.0\t1.0e-50\t!\trRNA
    chr2\t-\tmir-1\tRF004\tcm\t1\t70\t500\t570\t+\t-\t1\t0.43\t0.0\t55.0\t1.0e-1\t!\tmiRNA low
""")


@pytest.fixture
def tblout_file(tmp_path):
    f = tmp_path / 'search.tblout'
    f.write_text(TBLOUT_FILE)
    return str(f)


class TestParseTblout:
    def test_all_pass_at_permissive_threshold(self, tblout_file):
        hits = parse_tblout(tblout_file, max_evalue=1.0, min_score=0.0)
        assert len(hits) == 4

    def test_evalue_filter(self, tblout_file):
        hits = parse_tblout(tblout_file, max_evalue=0.01, min_score=0.0)
        assert len(hits) == 3  # 1e-1 fails

    def test_score_filter(self, tblout_file):
        hits = parse_tblout(tblout_file, max_evalue=1.0, min_score=80.0)
        assert len(hits) == 2  # mir-21 (87.3) and LSU_rRNA (1750.0)

    def test_biotypes_assigned(self, tblout_file):
        hits = parse_tblout(tblout_file, max_evalue=1.0, min_score=0.0)
        bts = {h.cm_name: h.biotype for h in hits}
        assert bts['mir-21'] == 'miRNA'
        assert bts['SNORD14'] == 'snoRNA'
        assert bts['LSU_rRNA'] == 'rRNA'


# ---------------------------------------------------------------------------
# parse_blast_tabular
# ---------------------------------------------------------------------------

BLAST_TABLE = textwrap.dedent("""\
    hsa-mir-21\tchr1\t95.0\t88\t4\t0\t1\t88\t1000\t1088\t1.0e-20\t150.0
    hsa-mir-155\tchr2\t70.0\t70\t0\t0\t1\t70\t5000\t5070\t1.0e-5\t80.0
    hsa-mir-200\tchr3\t50.0\t50\t0\t0\t1\t50\t100\t150\t0.5\t40.0
""")


@pytest.fixture
def blast_file(tmp_path):
    f = tmp_path / 'mirna.blast'
    f.write_text(BLAST_TABLE)
    return str(f)


class TestParseBlastTabular:
    def test_default_thresholds(self, blast_file):
        hits = parse_blast_tabular(blast_file, min_pident=80.0, max_evalue=0.01)
        assert len(hits) == 1  # only mir-21 passes both

    def test_permissive_returns_more(self, blast_file):
        hits = parse_blast_tabular(blast_file, min_pident=50.0, max_evalue=1.0)
        assert len(hits) == 3

    def test_pident_filter_only(self, blast_file):
        hits = parse_blast_tabular(blast_file, min_pident=80.0, max_evalue=1.0)
        assert len(hits) == 1


# ---------------------------------------------------------------------------
# write_gff3
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def _make_rfam_hit(self, seqname='chr1', cm_name='mir-21', seq_from=1000,
                       seq_to=1088, strand='+', score=87.3, evalue=1e-16):
        h = CmsearchHit(seqname=seqname, cm_name=cm_name, cm_acc='RF001',
                        seq_from=seq_from, seq_to=seq_to, strand=strand,
                        score=score, evalue=evalue)
        h.biotype = rfam_biotype(cm_name)
        return h

    def _make_blast_hit(self):
        return BlastHit('hsa-mir-21', 'chr2', 95.0, 5000, 5088, 1e-20, 150.0, '+')

    def test_rfam_hit_written(self, tmp_path):
        out = str(tmp_path / 'out.gff3')
        counts = write_gff3([self._make_rfam_hit()], [], out, 'test')
        assert counts.get('miRNA', 0) == 1
        content = open(out).read()
        assert '\tgene\t' in content
        assert '\ttranscript\t' in content
        assert '\texon\t' in content

    def test_blast_hit_written(self, tmp_path):
        out = str(tmp_path / 'out.gff3')
        counts = write_gff3([], [self._make_blast_hit()], out, 'test')
        assert counts.get('miRNA', 0) == 1

    def test_biotype_in_attrs(self, tmp_path):
        out = str(tmp_path / 'out.gff3')
        write_gff3([self._make_rfam_hit()], [], out, 'test')
        content = open(out).read()
        assert 'biotype=miRNA' in content

    def test_both_sources_combined(self, tmp_path):
        out = str(tmp_path / 'out.gff3')
        counts = write_gff3(
            [self._make_rfam_hit('chr1', 'SNORD14', 1000, 1065)],
            [self._make_blast_hit()],
            out, 'test'
        )
        assert sum(counts.values()) == 2

    def test_gff3_header(self, tmp_path):
        out = str(tmp_path / 'out.gff3')
        write_gff3([], [], out, 'test')
        assert open(out).readline().startswith('##gff-version 3')
