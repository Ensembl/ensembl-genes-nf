"""
Unit tests for load_gff3_to_core.py

Uses SQLite (via pymysql-compatible interface via a thin adapter) rather than
MySQL, since the schema DDL subset we need works identically on SQLite.
We use the standard `sqlite3` module directly in tests, and test the
Python functions in isolation.
"""

import json
import sys
import os
import textwrap
from pathlib import Path
from typing import Dict, List

import pytest

# Add bin/ to path
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))

from load_gff3_to_core import (
    GffExon, GffCds, GffTranscript, GffGene,
    parse_gff3,
    compute_translation_coords,
    exon_phase_from_cds,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def make_coding_transcript(
    tid="tx1", gid="gene1", seq="chr1",
    tx_start=100, tx_end=900, strand=1,
    exon_coords=((100, 300), (500, 900)),
    cds_coords=((150, 300), (500, 800)),
    biotype="protein_coding",
) -> GffTranscript:
    tx = GffTranscript(
        transcript_id=tid, gene_id=gid,
        seq_region_name=seq, start=tx_start, end=tx_end,
        strand=strand, biotype=biotype,
    )
    for s, e in exon_coords:
        tx.exons.append(GffExon(seq_region_name=seq, start=s, end=e, strand=strand))
    cds_phase = 0
    for s, e in cds_coords:
        tx.cds_features.append(GffCds(seq_region_name=seq, start=s, end=e,
                                       strand=strand, phase=cds_phase))
        cds_phase = (cds_phase + (e - s + 1)) % 3
    return tx


MINIMAL_GFF3 = textwrap.dedent("""\
    ##gff-version 3
    chr1\tensembl\tgene\t100\t900\t.\t+\t.\tID=gene1;Name=BRCA1;biotype=protein_coding
    chr1\tensembl\tmRNA\t100\t900\t.\t+\t.\tID=tx1;Parent=gene1;biotype=protein_coding
    chr1\tensembl\texon\t100\t300\t.\t+\t.\tID=exon1;Parent=tx1
    chr1\tensembl\texon\t500\t900\t.\t+\t.\tID=exon2;Parent=tx1
    chr1\tensembl\tCDS\t150\t300\t.\t+\t0\tID=cds1;Parent=tx1
    chr1\tensembl\tCDS\t500\t800\t.\t+\t0\tID=cds2;Parent=tx1
    chr1\tensembl\tgene\t1000\t1500\t.\t-\t.\tID=gene2;Name=TP53;biotype=protein_coding
    chr1\tensembl\tmRNA\t1000\t1500\t.\t-\t.\tID=tx2;Parent=gene2;biotype=protein_coding
    chr1\tensembl\texon\t1000\t1200\t.\t-\t.\tID=exon3;Parent=tx2
    chr1\tensembl\texon\t1350\t1500\t.\t-\t.\tID=exon4;Parent=tx2
    chr1\tensembl\tCDS\t1050\t1200\t.\t-\t0\tID=cds3;Parent=tx2
    chr1\tensembl\tCDS\t1350\t1480\t.\t-\t0\tID=cds4;Parent=tx2
""")

NONCODING_GFF3 = textwrap.dedent("""\
    ##gff-version 3
    chr1\tensembl\tgene\t100\t500\t.\t+\t.\tID=nc_gene;Name=SNRNA1;biotype=snRNA
    chr1\tensembl\ttranscript\t100\t500\t.\t+\t.\tID=nc_tx;Parent=nc_gene;biotype=snRNA
    chr1\tensembl\texon\t100\t500\t.\t+\t.\tID=nc_exon;Parent=nc_tx
""")


@pytest.fixture
def minimal_gff3_file(tmp_path):
    f = tmp_path / "minimal.gff3"
    f.write_text(MINIMAL_GFF3)
    return str(f)


@pytest.fixture
def noncoding_gff3_file(tmp_path):
    f = tmp_path / "noncoding.gff3"
    f.write_text(NONCODING_GFF3)
    return str(f)


# ---------------------------------------------------------------------------
# Tests: GFF3 parsing
# ---------------------------------------------------------------------------

class TestParseGff3:

    def test_returns_list_of_genes(self, minimal_gff3_file):
        genes = parse_gff3(minimal_gff3_file)
        assert len(genes) == 2

    def test_gene_coordinates(self, minimal_gff3_file):
        genes = {g.gene_id: g for g in parse_gff3(minimal_gff3_file)}
        g = genes["gene1"]
        assert g.start == 100
        assert g.end == 900
        assert g.strand == 1
        assert g.biotype == "protein_coding"

    def test_transcript_attached_to_gene(self, minimal_gff3_file):
        genes = {g.gene_id: g for g in parse_gff3(minimal_gff3_file)}
        assert len(genes["gene1"].transcripts) == 1
        assert genes["gene1"].transcripts[0].transcript_id == "tx1"

    def test_exons_parsed(self, minimal_gff3_file):
        genes = {g.gene_id: g for g in parse_gff3(minimal_gff3_file)}
        tx = genes["gene1"].transcripts[0]
        assert len(tx.exons) == 2
        assert tx.exons[0].start == 100
        assert tx.exons[0].end == 300

    def test_cds_features_parsed(self, minimal_gff3_file):
        genes = {g.gene_id: g for g in parse_gff3(minimal_gff3_file)}
        tx = genes["gene1"].transcripts[0]
        assert len(tx.cds_features) == 2
        assert tx.is_coding

    def test_minus_strand_gene(self, minimal_gff3_file):
        genes = {g.gene_id: g for g in parse_gff3(minimal_gff3_file)}
        g = genes["gene2"]
        assert g.strand == -1

    def test_noncoding_transcript_has_no_cds(self, noncoding_gff3_file):
        genes = parse_gff3(noncoding_gff3_file)
        tx = genes[0].transcripts[0]
        assert not tx.is_coding
        assert tx.biotype == "snRNA"


# ---------------------------------------------------------------------------
# Tests: compute_translation_coords
# ---------------------------------------------------------------------------

class TestComputeTranslationCoords:

    def test_plus_strand_simple(self):
        """Simple 2-exon + strand transcript."""
        tx = make_coding_transcript(
            exon_coords=((100, 300), (500, 900)),
            cds_coords=((150, 300), (500, 800)),
        )
        result = compute_translation_coords(tx)
        assert result is not None
        start_exon, seq_start, end_exon, seq_end = result
        # CDS starts at 150 in exon 100-300
        assert start_exon.start == 100
        assert seq_start == 150 - 100 + 1  # = 51
        # CDS ends at 800 in exon 500-900
        assert end_exon.start == 500
        assert seq_end == 800 - 500 + 1  # = 301

    def test_plus_strand_cds_equals_exon(self):
        """CDS spans full exon."""
        tx = make_coding_transcript(
            exon_coords=((100, 200), (300, 400)),
            cds_coords=((100, 200), (300, 350)),
        )
        result = compute_translation_coords(tx)
        assert result is not None
        start_exon, seq_start, _, _ = result
        assert seq_start == 1  # CDS at exon start

    def test_minus_strand(self):
        """Minus strand: exon 1350-1500, CDS 1350-1480."""
        tx = make_coding_transcript(
            strand=-1,
            exon_coords=((1000, 1200), (1350, 1500)),
            cds_coords=((1050, 1200), (1350, 1480)),
        )
        result = compute_translation_coords(tx)
        assert result is not None
        start_exon, seq_start, end_exon, seq_end = result
        # On minus strand, start_exon is the one with HIGHEST genomic coords
        # CDS "start" (first translated base) = 1480 in exon 1350-1500
        # seq_start = 1500 - 1480 + 1 = 21
        assert start_exon.start == 1350
        assert seq_start == 1500 - 1480 + 1

    def test_returns_none_when_no_cds(self):
        tx = GffTranscript("nc1", "g1", "chr1", 1, 100, 1, "snRNA")
        tx.exons.append(GffExon("chr1", 1, 100, 1))
        assert compute_translation_coords(tx) is None

    def test_returns_none_when_no_exons(self):
        tx = GffTranscript("t1", "g1", "chr1", 100, 300, 1, "protein_coding")
        tx.cds_features.append(GffCds("chr1", 150, 300, 1, 0))
        # No exons → can't find start/end exon
        assert compute_translation_coords(tx) is None


# ---------------------------------------------------------------------------
# Tests: exon_phase_from_cds
# ---------------------------------------------------------------------------

class TestExonPhaseFromCds:

    def test_noncoding_exons_have_phase_minus1(self):
        tx = GffTranscript("nc", "g", "chr1", 1, 100, 1, "snRNA")
        tx.exons.append(GffExon("chr1", 1, 100, 1))
        phases = exon_phase_from_cds(tx)
        assert phases[(1, 100)] == (-1, -1)

    def test_single_exon_in_frame(self):
        """Single exon, CDS length = 9 (3 codons, ends at phase 0)."""
        tx = make_coding_transcript(
            exon_coords=((100, 200),),
            cds_coords=((100, 108),),  # 9 bp CDS
        )
        phases = exon_phase_from_cds(tx)
        # phase at start of CDS = 0; end_phase = 9 % 3 = 0
        assert phases[(100, 200)] == (0, 0)

    def test_two_exon_cds_phase_carried_over(self):
        """First CDS exon is 7 bp (phase 0, end_phase 1); second starts at phase 1."""
        tx = make_coding_transcript(
            exon_coords=((100, 200), (300, 500)),
            cds_coords=((100, 106), (300, 450)),  # 7 bp + 151 bp
        )
        phases = exon_phase_from_cds(tx)
        phase_e1 = phases[(100, 200)]
        phase_e2 = phases[(300, 500)]
        assert phase_e1[0] == 0        # starts in phase 0
        assert phase_e1[1] == 7 % 3   # end_phase = 1
        assert phase_e2[0] == 7 % 3   # picks up from end of exon 1

    def test_fully_utr_exons_are_minus1(self):
        """First and last exons are purely UTR."""
        tx = make_coding_transcript(
            exon_coords=((50, 99), (100, 200), (300, 500), (501, 600)),
            cds_coords=((100, 200), (300, 450)),
        )
        phases = exon_phase_from_cds(tx)
        assert phases[(50, 99)]   == (-1, -1)   # 5' UTR exon
        assert phases[(501, 600)] == (-1, -1)   # 3' UTR exon
        assert phases[(100, 200)][0] == 0       # first coding exon phase = 0


# ---------------------------------------------------------------------------
# Tests: stable ID format
# ---------------------------------------------------------------------------

class TestStableIdFormat:
    """Test the stable ID string format (no DB needed)."""

    def _make_id(self, prefix, id_type, n):
        return f"ENS{prefix}{id_type}{str(n).zfill(11)}"

    def test_human_gene_id(self):
        assert self._make_id("", "G", 1) == "ENSG00000000001"

    def test_chicken_gene_id(self):
        assert self._make_id("GAL", "G", 42) == "ENSGALG00000000042"

    def test_transcript_id(self):
        assert self._make_id("", "T", 100) == "ENST00000000100"

    def test_exon_id(self):
        assert self._make_id("", "E", 999) == "ENSE00000000999"

    def test_translation_id(self):
        assert self._make_id("", "P", 12345678901) == "ENSP12345678901"
