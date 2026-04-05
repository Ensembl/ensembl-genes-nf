"""
Tests for add_utrs.py

Covers:
  - parse_gff3:       minimal GFF3 parsing (gene/mRNA/exon/CDS)
  - cds_matches_donor:  matching logic for various CDS/exon configurations
  - add_utr_to_transcript: UTR grafting, clipping, and filtering
  - write_gff3:       output structure
  - Integration:      end-to-end parse → match → graft → write
"""

from __future__ import annotations

import os
import sys
import textwrap

import pytest

# ---------------------------------------------------------------------------
# Path setup — allow importing from the sibling bin/ directory
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from add_utrs import (
    CdsExon,
    Exon,
    Gene,
    Transcript,
    _build_donor_index,
    _clip_utr,
    _parse_attrs,
    _set_attr,
    add_utr_to_transcript,
    cds_matches_donor,
    find_matching_donor,
    parse_gff3,
    write_gff3,
)


# ===========================================================================
# Helper factories
# ===========================================================================

def make_transcript(
    seqname: str = 'chr1',
    start: int = 1000,
    end: int = 5000,
    strand: str = '+',
    tx_id: str = 'tx1',
    gene_id: str = 'gene1',
    exons=None,   # list of (start, end)
    cds=None,     # list of (start, end) — phase defaults to '0'
    attrs: str = '',
) -> Transcript:
    tx = Transcript(
        seqname=seqname,
        start=start,
        end=end,
        strand=strand,
        tx_id=tx_id,
        gene_id=gene_id,
        source='test',
        score='.',
        attrs=attrs or f'ID={tx_id};Parent={gene_id}',
    )
    tx.exons = [Exon(start=s, end=e) for s, e in (exons or [])]
    tx.cds   = [CdsExon(start=s, end=e, phase='0') for s, e in (cds or [])]
    return tx


def make_gene(
    gene_id: str = 'gene1',
    seqname: str = 'chr1',
    start: int = 1000,
    end: int = 5000,
    strand: str = '+',
    transcripts=None,
) -> Gene:
    g = Gene(
        seqname=seqname,
        start=start,
        end=end,
        strand=strand,
        gene_id=gene_id,
        source='test',
        score='.',
        attrs=f'ID={gene_id}',
    )
    g.transcripts = transcripts or []
    return g


def write_tmp_gff3(tmp_path, filename: str, content: str) -> str:
    p = tmp_path / filename
    p.write_text(textwrap.dedent(content))
    return str(p)


# ===========================================================================
# _parse_attrs
# ===========================================================================

class TestParseAttrs:
    def test_simple_key_value(self):
        d = _parse_attrs('ID=gene1;Parent=chr1')
        assert d['ID'] == 'gene1'
        assert d['Parent'] == 'chr1'

    def test_single_key(self):
        d = _parse_attrs('ID=tx1')
        assert d == {'ID': 'tx1'}

    def test_empty_string(self):
        d = _parse_attrs('')
        assert d == {}

    def test_trailing_semicolon(self):
        d = _parse_attrs('ID=tx1;biotype=protein_coding;')
        assert d['ID'] == 'tx1'
        assert d['biotype'] == 'protein_coding'


# ===========================================================================
# _set_attr
# ===========================================================================

class TestSetAttr:
    def test_appends_new_key(self):
        result = _set_attr('ID=tx1', 'utr_source', 'donor.gff3')
        assert 'utr_source=donor.gff3' in result
        assert 'ID=tx1' in result

    def test_replaces_existing_key(self):
        result = _set_attr('ID=tx1;utr_source=old.gff3', 'utr_source', 'new.gff3')
        assert 'utr_source=new.gff3' in result
        assert 'utr_source=old.gff3' not in result

    def test_other_keys_preserved(self):
        result = _set_attr('ID=tx1;biotype=protein_coding', 'utr_source', 'donor.gff3')
        assert 'ID=tx1' in result
        assert 'biotype=protein_coding' in result


# ===========================================================================
# parse_gff3
# ===========================================================================

MINIMAL_GFF3 = """\
    ##gff-version 3
    chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=gene1;biotype=protein_coding
    chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding
    chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tParent=tx1
    chr1\tsrc\texon\t3000\t5000\t.\t+\t.\tParent=tx1
    chr1\tsrc\tCDS\t1050\t2000\t.\t+\t0\tParent=tx1
    chr1\tsrc\tCDS\t3000\t4950\t.\t+\t0\tParent=tx1
    """

TWO_GENE_GFF3 = """\
    ##gff-version 3
    chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=gene1;biotype=protein_coding
    chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding
    chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tParent=tx1
    chr1\tsrc\texon\t3000\t5000\t.\t+\t.\tParent=tx1
    chr1\tsrc\tCDS\t1050\t2000\t.\t+\t0\tParent=tx1
    chr1\tsrc\tCDS\t3000\t4950\t.\t+\t0\tParent=tx1
    chr2\tsrc\tgene\t100\t900\t.\t-\t.\tID=gene2;biotype=protein_coding
    chr2\tsrc\tmRNA\t100\t900\t.\t-\t.\tID=tx2;Parent=gene2;biotype=protein_coding
    chr2\tsrc\texon\t100\t400\t.\t-\t.\tParent=tx2
    chr2\tsrc\texon\t600\t900\t.\t-\t.\tParent=tx2
    chr2\tsrc\tCDS\t150\t400\t.\t-\t0\tParent=tx2
    chr2\tsrc\tCDS\t600\t850\t.\t-\t0\tParent=tx2
    """

DONOR_GFF3 = """\
    ##gff-version 3
    chr1\tdonor\tgene\t800\t5200\t.\t+\t.\tID=dgene1
    chr1\tdonor\tmRNA\t800\t5200\t.\t+\t.\tID=dtx1;Parent=dgene1
    chr1\tdonor\texon\t800\t2000\t.\t+\t.\tParent=dtx1
    chr1\tdonor\texon\t3000\t5200\t.\t+\t.\tParent=dtx1
    """


class TestParseGff3:
    def test_parses_gene_and_transcript(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        assert 'gene1' in genes
        assert len(genes['gene1'].transcripts) == 1

    def test_transcript_has_exons(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        tx = genes['gene1'].transcripts[0]
        assert len(tx.exons) == 2

    def test_transcript_has_cds(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        tx = genes['gene1'].transcripts[0]
        assert len(tx.cds) == 2

    def test_exon_coordinates(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        tx = genes['gene1'].transcripts[0]
        assert tx.exons[0].start == 1000
        assert tx.exons[0].end   == 2000
        assert tx.exons[1].start == 3000
        assert tx.exons[1].end   == 5000

    def test_cds_coordinates(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        tx = genes['gene1'].transcripts[0]
        assert tx.cds[0].start == 1050
        assert tx.cds[1].end   == 4950

    def test_strand_preserved(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'minimal.gff3', MINIMAL_GFF3)
        genes = parse_gff3(path)
        assert genes['gene1'].transcripts[0].strand == '+'

    def test_two_genes(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'two.gff3', TWO_GENE_GFF3)
        genes = parse_gff3(path)
        assert len(genes) == 2
        assert 'gene1' in genes
        assert 'gene2' in genes

    def test_minus_strand(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'two.gff3', TWO_GENE_GFF3)
        genes = parse_gff3(path)
        assert genes['gene2'].transcripts[0].strand == '-'

    def test_empty_file(self, tmp_path):
        path = write_tmp_gff3(tmp_path, 'empty.gff3', '##gff-version 3\n')
        genes = parse_gff3(path)
        assert genes == {}

    def test_exons_sorted(self, tmp_path):
        # Write exons in reverse order to test sorting
        content = """\
            ##gff-version 3
            chr1\tsrc\tgene\t1\t5000\t.\t+\t.\tID=g1
            chr1\tsrc\tmRNA\t1\t5000\t.\t+\t.\tID=t1;Parent=g1
            chr1\tsrc\texon\t3000\t5000\t.\t+\t.\tParent=t1
            chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tParent=t1
            chr1\tsrc\tCDS\t1000\t2000\t.\t+\t0\tParent=t1
            """
        path = write_tmp_gff3(tmp_path, 'unsorted.gff3', content)
        genes = parse_gff3(path)
        tx = genes['g1'].transcripts[0]
        starts = [e.start for e in tx.exons]
        assert starts == sorted(starts)


# ===========================================================================
# cds_matches_donor
# ===========================================================================

class TestCdsMatchesDonor:
    def _cds(self, *intervals):
        return [CdsExon(start=s, end=e, phase='0') for s, e in intervals]

    def _exons(self, *intervals):
        return [Exon(start=s, end=e) for s, e in intervals]

    # -- Positive cases ------------------------------------------------------

    def test_exact_cds_match(self):
        """Donor exons exactly equal CDS exons."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((1000, 2000), (3000, 5000))
        assert cds_matches_donor(cds, donor) is True

    def test_five_prime_utr_extension(self):
        """Donor extends the first exon to the left (5' UTR)."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((800, 2000), (3000, 5000))
        assert cds_matches_donor(cds, donor) is True

    def test_three_prime_utr_extension(self):
        """Donor extends the last exon to the right (3' UTR)."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((1000, 2000), (3000, 5200))
        assert cds_matches_donor(cds, donor) is True

    def test_both_utr_extensions(self):
        """Donor extends both the first and last exon."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((800, 2000), (3000, 5200))
        assert cds_matches_donor(cds, donor) is True

    def test_utr_as_separate_exon_5prime(self):
        """Donor has a dedicated 5' UTR exon separate from the CDS exon."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((700, 900), (1000, 2000), (3000, 5000))
        assert cds_matches_donor(cds, donor) is True

    def test_utr_as_separate_exon_3prime(self):
        """Donor has a dedicated 3' UTR exon separate from the last CDS exon."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((1000, 2000), (3000, 5000), (5500, 5800))
        assert cds_matches_donor(cds, donor) is True

    def test_single_exon_cds_exact(self):
        """Single-exon CDS matched by a single donor exon."""
        cds   = self._cds((1000, 5000))
        donor = self._exons((1000, 5000))
        assert cds_matches_donor(cds, donor) is True

    def test_single_exon_cds_with_both_utrs(self):
        """Donor spans more than the single-exon CDS."""
        cds   = self._cds((1000, 5000))
        donor = self._exons((800, 5200))
        assert cds_matches_donor(cds, donor) is True

    def test_three_exon_cds_all_junctions_match(self):
        """Three-exon CDS: both internal junctions must match."""
        cds   = self._cds((1000, 2000), (3000, 4000), (5000, 6000))
        donor = self._exons((800, 2000), (3000, 4000), (5000, 6200))
        assert cds_matches_donor(cds, donor) is True

    # -- Negative cases ------------------------------------------------------

    def test_different_internal_splice_site(self):
        """Donor has a different internal exon boundary — no match."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        # Donor ends at 2100 instead of 2000 — internal boundary differs
        donor = self._exons((1000, 2100), (3000, 5000))
        assert cds_matches_donor(cds, donor) is False

    def test_different_internal_start(self):
        """Donor's second exon starts at wrong position."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((1000, 2000), (3100, 5000))
        assert cds_matches_donor(cds, donor) is False

    def test_donor_does_not_cover_cds_start(self):
        """Donor exon starts after the CDS start — does not contain CDS."""
        cds   = self._cds((1000, 2000), (3000, 5000))
        donor = self._exons((1100, 2000), (3000, 5000))
        assert cds_matches_donor(cds, donor) is False

    def test_empty_cds(self):
        """Empty CDS list returns False."""
        donor = self._exons((1000, 5000))
        assert cds_matches_donor([], donor) is False

    def test_empty_donor(self):
        """Empty donor list returns False."""
        cds = self._cds((1000, 5000))
        assert cds_matches_donor(cds, []) is False

    def test_single_exon_donor_does_not_span_cds(self):
        """Single-exon donor that doesn't fully span the single-exon CDS."""
        cds   = self._cds((1000, 5000))
        donor = self._exons((1000, 4500))
        assert cds_matches_donor(cds, donor) is False


# ===========================================================================
# _clip_utr
# ===========================================================================

class TestClipUtr:
    def _exons(self, *intervals):
        return [Exon(start=s, end=e) for s, e in intervals]

    def test_no_clipping_needed(self):
        """Total length within max — no change."""
        exons = self._exons((100, 199))   # 100 bp
        result = _clip_utr(exons, max_len=200, from_right=False)
        assert result == exons

    def test_clip_3prime_single_exon(self):
        """Single 3' UTR exon clipped to max_len."""
        exons = self._exons((5001, 6000))   # 1000 bp
        result = _clip_utr(exons, max_len=500, from_right=False)
        assert len(result) == 1
        assert result[0].end - result[0].start + 1 == 500

    def test_clip_5prime_single_exon(self):
        """Single 5' UTR exon clipped from the left (from_right=True)."""
        exons = self._exons((1, 1000))   # 1000 bp, adjacent to CDS at 1001
        result = _clip_utr(exons, max_len=200, from_right=True)
        assert len(result) == 1
        assert result[0].end - result[0].start + 1 == 200
        # CDS-proximal end (right side) must be preserved
        assert result[0].end == 1000

    def test_clip_5prime_two_exons_removes_distal(self):
        """
        5' UTR has two exons; max allows only the proximal one.
        exons (sorted ascending): (500, 699) [200 bp], (800, 999) [200 bp]
        CDS starts at 1000. from_right=True means we iterate from right.
        """
        exons = self._exons((500, 699), (800, 999))   # 200 + 200 = 400 bp
        result = _clip_utr(exons, max_len=200, from_right=True)
        # Should keep only the proximal exon (800-999)
        assert len(result) == 1
        assert result[0].start == 800

    def test_empty_utr(self):
        result = _clip_utr([], max_len=500, from_right=False)
        assert result == []


# ===========================================================================
# add_utr_to_transcript
# ===========================================================================

class TestAddUtrToTranscript:
    """Tests for the UTR grafting function."""

    def _make_acceptor(self, cds_intervals):
        """Minimal acceptor with given CDS, no pre-existing exons."""
        starts = [s for s, _ in cds_intervals]
        ends   = [e for _, e in cds_intervals]
        tx = make_transcript(
            start=min(starts),
            end=max(ends),
            cds=cds_intervals,
            exons=cds_intervals,  # acceptor exons = CDS only (no UTR)
        )
        return tx

    def _make_donor(self, exon_intervals):
        """Donor transcript with given exon structure."""
        starts = [s for s, _ in exon_intervals]
        ends   = [e for _, e in exon_intervals]
        return make_transcript(
            start=min(starts),
            end=max(ends),
            exons=exon_intervals,
        )

    def test_adds_5prime_utr_exon(self):
        """Donor has a dedicated 5' UTR exon — it should appear in output."""
        acceptor = self._make_acceptor([(1000, 2000), (3000, 5000)])
        donor    = self._make_donor([(700, 999), (1000, 2000), (3000, 5000)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=1)
        starts = [e.start for e in result]
        assert 700 in starts

    def test_adds_3prime_utr_exon(self):
        """Donor has a dedicated 3' UTR exon — it should appear in output."""
        acceptor = self._make_acceptor([(1000, 2000), (3000, 5000)])
        donor    = self._make_donor([(1000, 2000), (3000, 5000), (5200, 5500)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=1)
        ends = [e.end for e in result]
        assert 5500 in ends

    def test_cds_exon_boundaries_preserved(self):
        """CDS exon coordinates must not be changed by UTR grafting."""
        acceptor = self._make_acceptor([(1000, 2000), (3000, 5000)])
        donor    = self._make_donor([(700, 2000), (3000, 5200)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=1)
        # Find the exons covering the CDS
        exon_starts = {e.start for e in result}
        exon_ends   = {e.end   for e in result}
        assert 1000 in exon_starts
        assert 2000 in exon_ends
        assert 3000 in exon_starts
        assert 5000 in exon_ends

    def test_respects_max_5prime_limit(self):
        """5' UTR extension must be capped at max_5prime bp."""
        # Donor has a 600 bp 5' UTR exon; limit is 200 bp
        acceptor = self._make_acceptor([(1000, 2000)])
        donor    = self._make_donor([(400, 999), (1000, 2000)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=200, max_3prime=10000,
                                       min_exon=1)
        five_prime = [e for e in result if e.end < 1000]
        total_5utr = sum(e.end - e.start + 1 for e in five_prime)
        assert total_5utr <= 200

    def test_respects_max_3prime_limit(self):
        """3' UTR extension must be capped at max_3prime bp."""
        # Donor extends 1000 bp beyond CDS end; limit is 300 bp
        acceptor = self._make_acceptor([(1000, 2000)])
        donor    = self._make_donor([(1000, 2000), (2001, 3000)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=300,
                                       min_exon=1)
        three_prime = [e for e in result if e.start > 2000]
        total_3utr = sum(e.end - e.start + 1 for e in three_prime)
        assert total_3utr <= 300

    def test_min_utr_exon_size_filter(self):
        """UTR exons smaller than min_exon are discarded."""
        # Donor has a tiny 5' UTR exon of 10 bp; min_exon=30 → should be dropped
        acceptor = self._make_acceptor([(1000, 2000)])
        donor    = self._make_donor([(990, 999), (1000, 2000)])  # 10 bp UTR exon
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=30)
        five_prime = [e for e in result if e.end < 1000]
        assert five_prime == []

    def test_min_utr_exon_size_accepts_large_enough(self):
        """UTR exons >= min_exon are kept."""
        acceptor = self._make_acceptor([(1000, 2000)])
        donor    = self._make_donor([(900, 999), (1000, 2000)])  # 100 bp UTR exon
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=30)
        five_prime = [e for e in result if e.end < 1000]
        assert len(five_prime) == 1

    def test_no_donor_utr_returns_cds_exons(self):
        """When donor exons equal CDS exactly, output should still have CDS exons."""
        acceptor = self._make_acceptor([(1000, 2000), (3000, 5000)])
        donor    = self._make_donor([(1000, 2000), (3000, 5000)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=30)
        assert len(result) == 2

    def test_output_exons_sorted(self):
        """Output exons must be sorted by start position."""
        acceptor = self._make_acceptor([(1000, 2000), (3000, 5000)])
        donor    = self._make_donor([(700, 999), (1000, 2000), (3000, 5000), (5100, 5300)])
        result = add_utr_to_transcript(acceptor, donor,
                                       max_5prime=5000, max_3prime=10000,
                                       min_exon=1)
        starts = [e.start for e in result]
        assert starts == sorted(starts)


# ===========================================================================
# find_matching_donor
# ===========================================================================

class TestFindMatchingDonor:
    def test_finds_matching_donor(self):
        """Donor transcript that contains CDS must be found."""
        acceptor = make_transcript(
            seqname='chr1', strand='+',
            cds=[(1000, 2000), (3000, 5000)],
        )
        donor_tx = make_transcript(
            tx_id='dtx1', gene_id='dgene1',
            seqname='chr1', strand='+',
            exons=[(800, 2000), (3000, 5200)],
        )
        gene = make_gene(gene_id='dgene1', transcripts=[donor_tx])
        index = _build_donor_index({'dgene1': gene})
        result = find_matching_donor(acceptor, index)
        assert result is donor_tx

    def test_no_match_wrong_chromosome(self):
        """Donor on different chromosome must not match."""
        acceptor = make_transcript(seqname='chr1', strand='+',
                                   cds=[(1000, 2000), (3000, 5000)])
        donor_tx = make_transcript(tx_id='dtx1', gene_id='dg1',
                                   seqname='chr2', strand='+',
                                   exons=[(800, 2000), (3000, 5200)])
        gene = make_gene(gene_id='dg1', seqname='chr2', transcripts=[donor_tx])
        index = _build_donor_index({'dg1': gene})
        result = find_matching_donor(acceptor, index)
        assert result is None

    def test_no_match_wrong_strand(self):
        """Donor on opposite strand must not match."""
        acceptor = make_transcript(seqname='chr1', strand='+',
                                   cds=[(1000, 2000), (3000, 5000)])
        donor_tx = make_transcript(tx_id='dtx1', gene_id='dg1',
                                   seqname='chr1', strand='-',
                                   exons=[(800, 2000), (3000, 5200)])
        gene = make_gene(gene_id='dg1', transcripts=[donor_tx])
        index = _build_donor_index({'dg1': gene})
        result = find_matching_donor(acceptor, index)
        assert result is None

    def test_no_match_different_junctions(self):
        """Donor with incompatible splice sites must not match."""
        acceptor = make_transcript(seqname='chr1', strand='+',
                                   cds=[(1000, 2000), (3000, 5000)])
        donor_tx = make_transcript(tx_id='dtx1', gene_id='dg1',
                                   seqname='chr1', strand='+',
                                   exons=[(800, 2100), (3000, 5200)])
        gene = make_gene(gene_id='dg1', transcripts=[donor_tx])
        index = _build_donor_index({'dg1': gene})
        result = find_matching_donor(acceptor, index)
        assert result is None

    def test_acceptor_without_cds_returns_none(self):
        """Acceptor with no CDS must return None immediately."""
        acceptor = make_transcript(seqname='chr1', strand='+', cds=[])
        index = {}
        result = find_matching_donor(acceptor, index)
        assert result is None


# ===========================================================================
# write_gff3
# ===========================================================================

class TestWriteGff3:
    def test_writes_gff_header(self, tmp_path):
        gene = make_gene()
        tx = make_transcript(exons=[(1000, 2000)], cds=[(1000, 2000)])
        gene.transcripts = [tx]
        out = str(tmp_path / 'out.gff3')
        write_gff3({'gene1': gene}, out)
        with open(out) as fh:
            first = fh.readline().strip()
        assert first == '##gff-version 3'

    def test_writes_gene_line(self, tmp_path):
        gene = make_gene()
        tx = make_transcript(exons=[(1000, 2000)], cds=[(1000, 2000)])
        gene.transcripts = [tx]
        out = str(tmp_path / 'out.gff3')
        write_gff3({'gene1': gene}, out)
        content = open(out).read()
        assert '\tgene\t' in content

    def test_writes_mrna_line(self, tmp_path):
        gene = make_gene()
        tx = make_transcript(exons=[(1000, 2000)], cds=[(1000, 2000)])
        gene.transcripts = [tx]
        out = str(tmp_path / 'out.gff3')
        write_gff3({'gene1': gene}, out)
        content = open(out).read()
        assert '\tmRNA\t' in content

    def test_writes_exon_lines(self, tmp_path):
        gene = make_gene()
        tx = make_transcript(exons=[(1000, 2000), (3000, 4000)],
                             cds=[(1000, 2000), (3000, 4000)])
        gene.transcripts = [tx]
        out = str(tmp_path / 'out.gff3')
        write_gff3({'gene1': gene}, out)
        content = open(out).read()
        exon_lines = [l for l in content.splitlines() if '\texon\t' in l]
        assert len(exon_lines) == 2

    def test_writes_cds_lines(self, tmp_path):
        gene = make_gene()
        tx = make_transcript(exons=[(1000, 2000)], cds=[(1100, 1900)])
        gene.transcripts = [tx]
        out = str(tmp_path / 'out.gff3')
        write_gff3({'gene1': gene}, out)
        content = open(out).read()
        assert '\tCDS\t' in content

    def test_returns_correct_counts(self, tmp_path):
        gene = make_gene()
        tx1 = make_transcript(tx_id='tx1', exons=[(1000, 2000)], cds=[(1000, 2000)])
        tx2 = make_transcript(tx_id='tx2', start=3000, end=5000,
                              exons=[(3000, 5000)], cds=[(3000, 5000)])
        gene.transcripts = [tx1, tx2]
        out = str(tmp_path / 'out.gff3')
        g_count, tx_count = write_gff3({'gene1': gene}, out)
        assert g_count == 1
        assert tx_count == 2


# ===========================================================================
# Integration test
# ===========================================================================

CONSOLIDATED_GFF3 = """\
##gff-version 3
chr1\tconsolidation\tgene\t1000\t5000\t.\t+\t.\tID=cons_gene1;biotype=protein_coding
chr1\tconsolidation\tmRNA\t1000\t5000\t.\t+\t.\tID=cons_tx1;Parent=cons_gene1;biotype=protein_coding
chr1\tconsolidation\texon\t1000\t2000\t.\t+\t.\tParent=cons_tx1
chr1\tconsolidation\texon\t3000\t5000\t.\t+\t.\tParent=cons_tx1
chr1\tconsolidation\tCDS\t1050\t2000\t.\t+\t0\tParent=cons_tx1
chr1\tconsolidation\tCDS\t3000\t4950\t.\t+\t0\tParent=cons_tx1
chr2\tconsolidation\tgene\t100\t900\t.\t-\t.\tID=cons_gene2;biotype=protein_coding
chr2\tconsolidation\tmRNA\t100\t900\t.\t-\t.\tID=cons_tx2;Parent=cons_gene2;biotype=protein_coding
chr2\tconsolidation\texon\t150\t400\t.\t-\t.\tParent=cons_tx2
chr2\tconsolidation\texon\t600\t850\t.\t-\t.\tParent=cons_tx2
chr2\tconsolidation\tCDS\t150\t400\t.\t-\t0\tParent=cons_tx2
chr2\tconsolidation\tCDS\t600\t850\t.\t-\t0\tParent=cons_tx2
"""

DONOR_WITH_UTR_GFF3 = """\
##gff-version 3
chr1\trnaseq\tgene\t800\t5200\t.\t+\t.\tID=donor_gene1
chr1\trnaseq\tmRNA\t800\t5200\t.\t+\t.\tID=donor_tx1;Parent=donor_gene1
chr1\trnaseq\texon\t800\t2000\t.\t+\t.\tParent=donor_tx1
chr1\trnaseq\texon\t3000\t5200\t.\t+\t.\tParent=donor_tx1
"""


class TestIntegration:
    def test_utr_added_to_matched_transcript(self, tmp_path):
        """
        Full pipeline: parse acceptors + donor, match, graft UTRs, write output.
        The chr1 acceptor should receive UTR from the donor.
        """
        consolidated_path = write_tmp_gff3(
            tmp_path, 'cons.gff3', CONSOLIDATED_GFF3
        )
        donor_path = write_tmp_gff3(
            tmp_path, 'donor.gff3', DONOR_WITH_UTR_GFF3
        )
        out_path = str(tmp_path / 'out.with_utrs.gff3')

        # Parse
        acceptor_genes = parse_gff3(consolidated_path)
        donor_genes    = parse_gff3(donor_path)
        donor_index    = _build_donor_index(donor_genes)

        # Match and graft
        for gene in acceptor_genes.values():
            for tx in gene.transcripts:
                if not tx.cds:
                    continue
                matched = find_matching_donor(tx, donor_index)
                if matched:
                    new_exons = add_utr_to_transcript(
                        tx, matched,
                        max_5prime=5000, max_3prime=10000, min_exon=10
                    )
                    tx.exons  = new_exons
                    tx.start  = new_exons[0].start
                    tx.end    = new_exons[-1].end
                    tx.attrs  = _set_attr(tx.attrs, 'utr_source', 'donor.gff3')

        write_gff3(acceptor_genes, out_path)

        content = open(out_path).read()
        assert '##gff-version 3' in content

        # The chr1 transcript should now have a 5' UTR exon (starts before 1050)
        lines = content.splitlines()
        exon_lines = [l for l in lines if '\texon\t' in l]
        exon_starts = [int(l.split('\t')[3]) for l in exon_lines]
        # Donor 5' UTR exon starts at 800 — should appear
        assert any(s <= 900 for s in exon_starts), (
            f"Expected a 5' UTR exon starting at ~800; found starts: {sorted(exon_starts)}"
        )

    def test_chr2_unmatched_transcript_preserved(self, tmp_path):
        """
        chr2 transcript has no matching donor — it must still appear in output.
        """
        consolidated_path = write_tmp_gff3(
            tmp_path, 'cons.gff3', CONSOLIDATED_GFF3
        )
        donor_path = write_tmp_gff3(
            tmp_path, 'donor.gff3', DONOR_WITH_UTR_GFF3
        )
        out_path = str(tmp_path / 'out.with_utrs.gff3')

        acceptor_genes = parse_gff3(consolidated_path)
        donor_genes    = parse_gff3(donor_path)
        donor_index    = _build_donor_index(donor_genes)

        for gene in acceptor_genes.values():
            for tx in gene.transcripts:
                if not tx.cds:
                    continue
                matched = find_matching_donor(tx, donor_index)
                if matched:
                    new_exons = add_utr_to_transcript(
                        tx, matched,
                        max_5prime=5000, max_3prime=10000, min_exon=10
                    )
                    tx.exons = new_exons
                    tx.start = new_exons[0].start
                    tx.end   = new_exons[-1].end
                else:
                    # No match — ensure exons are populated from CDS
                    if not tx.exons:
                        tx.exons = [Exon(start=c.start, end=c.end) for c in tx.cds]

        write_gff3(acceptor_genes, out_path)
        content = open(out_path).read()
        assert 'cons_gene2' in content or 'cons_tx2' in content

    def test_utr_source_attribute_set(self, tmp_path):
        """Modified transcript must carry utr_source= attribute."""
        consolidated_path = write_tmp_gff3(
            tmp_path, 'cons.gff3', CONSOLIDATED_GFF3
        )
        donor_path = write_tmp_gff3(
            tmp_path, 'donor.gff3', DONOR_WITH_UTR_GFF3
        )
        out_path = str(tmp_path / 'out.with_utrs.gff3')

        acceptor_genes = parse_gff3(consolidated_path)
        donor_genes    = parse_gff3(donor_path)
        donor_index    = _build_donor_index(donor_genes)

        for gene in acceptor_genes.values():
            for tx in gene.transcripts:
                if not tx.cds:
                    continue
                matched = find_matching_donor(tx, donor_index)
                if matched:
                    new_exons = add_utr_to_transcript(
                        tx, matched,
                        max_5prime=5000, max_3prime=10000, min_exon=10
                    )
                    tx.exons = new_exons
                    tx.start = new_exons[0].start
                    tx.end   = new_exons[-1].end
                    tx.attrs = _set_attr(tx.attrs, 'utr_source', 'donor.gff3')

        write_gff3(acceptor_genes, out_path)
        content = open(out_path).read()
        assert 'utr_source=donor.gff3' in content
