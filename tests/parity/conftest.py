"""
Shared fixtures for biological parity tests.

Provides synthetic GFF3 files that model different annotation scenarios:
  - reference_gff3: the "ground truth" (mimics Perl/eHive output)
  - test_gff3: the "new pipeline" output, identical to reference but with
    controlled differences introduced to verify detection of regressions/improvements

Scenario families:
  1. identical  — both pipelines agree exactly
  2. utr_added  — new pipeline adds 5'/3' UTR that reference lacks
  3. cds_shift  — new pipeline has a shifted CDS boundary (regression)
  4. novel_gene — new pipeline calls an extra gene not in reference
  5. missed_gene — new pipeline misses a gene present in reference
  6. pseudo_recall — pseudogene flagged in reference but not in new pipeline
"""

import textwrap
from pathlib import Path

import pytest

PARITY_FIXTURES = Path(__file__).parent / "fixtures"


def _write_gff3(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content))
    return path


# ---------------------------------------------------------------------------
# Helpers for building GFF3 content strings
# ---------------------------------------------------------------------------

def _gene_block(
    chrom="chr1", gene_start=1000, gene_end=2000,
    strand="+", gene_id="gene1", biotype="protein_coding",
    tx_id="tx1",
    exons=((1000, 1200), (1600, 2000)),
    cds=((1050, 1200), (1600, 1950)),
    utr5=None,      # list of (start, end) or None
    utr3=None,
    is_canonical=True,
    is_pseudogene=False,
) -> str:
    feature_type = "mRNA" if not is_pseudogene else "pseudogenic_transcript"
    biotype_attr = f";biotype={biotype}" if biotype else ""
    canonical_attr = ";canonical_transcript=1" if is_canonical else ""

    lines = [
        f"{chrom}\t.\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={gene_id}{biotype_attr}",
        f"{chrom}\t.\t{feature_type}\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={tx_id};Parent={gene_id}{biotype_attr}{canonical_attr}",
    ]
    for i, (es, ee) in enumerate(exons, 1):
        lines.append(f"{chrom}\t.\texon\t{es}\t{ee}\t.\t{strand}\t.\tID={tx_id}.exon{i};Parent={tx_id}")
    for cs, ce in (cds or []):
        lines.append(f"{chrom}\t.\tCDS\t{cs}\t{ce}\t.\t{strand}\t0\tID={tx_id}.CDS;Parent={tx_id}")
    for us, ue in (utr5 or []):
        lines.append(f"{chrom}\t.\tfive_prime_UTR\t{us}\t{ue}\t.\t{strand}\t.\tID={tx_id}.utr5;Parent={tx_id}")
    for us, ue in (utr3 or []):
        lines.append(f"{chrom}\t.\tthree_prime_UTR\t{us}\t{ue}\t.\t{strand}\t.\tID={tx_id}.utr3;Parent={tx_id}")
    return "\n".join(lines)


def _header() -> str:
    return "##gff-version 3"


# ---------------------------------------------------------------------------
# Fixtures: identical annotation
# ---------------------------------------------------------------------------

@pytest.fixture()
def identical_ref_gff3(tmp_path) -> Path:
    """Reference GFF3 with two protein-coding genes, no UTRs."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="gene1", tx_id="tx1",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
        _gene_block(
            chrom="chr1", gene_start=3000, gene_end=4000, strand="-",
            gene_id="gene2", tx_id="tx2",
            exons=((3000, 3400), (3700, 4000)),
            cds=((3050, 3400), (3700, 3950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "ref_identical.gff3", content)


@pytest.fixture()
def identical_test_gff3(tmp_path) -> Path:
    """Test GFF3 identical to ref but with different gene IDs (as expected from new pipeline)."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="ENSG00000000001", tx_id="ENST00000000001",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
        _gene_block(
            chrom="chr1", gene_start=3000, gene_end=4000, strand="-",
            gene_id="ENSG00000000002", tx_id="ENST00000000002",
            exons=((3000, 3400), (3700, 4000)),
            cds=((3050, 3400), (3700, 3950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "test_identical.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: UTR added by new pipeline
# ---------------------------------------------------------------------------

@pytest.fixture()
def utr_ref_gff3(tmp_path) -> Path:
    """Reference with protein-coding gene, no UTR features."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="gene1", tx_id="tx1",
            exons=((1050, 1200), (1600, 1950)),
            cds=((1050, 1200), (1600, 1950)),
            utr5=None, utr3=None,
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "ref_utr.gff3", content)


@pytest.fixture()
def utr_test_gff3(tmp_path) -> Path:
    """Test GFF3 where utr_addition pipeline extended with 5' and 3' UTR."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="ENSG00000000001", tx_id="ENST00000000001",
            gene_start=900, gene_end=2100,
            exons=((900, 1200), (1600, 2100)),
            cds=((1050, 1200), (1600, 1950)),
            utr5=[(900, 1049)],   # 150 bp of 5' UTR
            utr3=[(1951, 2100)],  # 150 bp of 3' UTR
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "test_utr.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: CDS boundary shift (regression)
# ---------------------------------------------------------------------------

@pytest.fixture()
def shifted_ref_gff3(tmp_path) -> Path:
    """Reference gene with CDS at 1050-1200."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="gene1", tx_id="tx1",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "ref_shifted.gff3", content)


@pytest.fixture()
def shifted_test_gff3(tmp_path) -> Path:
    """Test gene with CDS shifted by 3 bp at start (exon 1 CDS starts at 1053 instead of 1050)."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="ENSG00000000001", tx_id="ENST00000000001",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1053, 1200), (1600, 1950)),   # shifted CDS start
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "test_shifted.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: novel gene in test (no reference match)
# ---------------------------------------------------------------------------

@pytest.fixture()
def novel_gene_ref_gff3(tmp_path) -> Path:
    """Reference with one gene."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="gene1", tx_id="tx1",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "ref_novel.gff3", content)


@pytest.fixture()
def novel_gene_test_gff3(tmp_path) -> Path:
    """Test GFF3 with the matched gene plus an entirely new gene on chr2."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="ENSG00000000001", tx_id="ENST00000000001",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
        _gene_block(
            chrom="chr2", gene_start=5000, gene_end=6000,
            gene_id="ENSG00000000002", tx_id="ENST00000000002",
            exons=((5000, 5300), (5700, 6000)),
            cds=((5050, 5300), (5700, 5950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "test_novel.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: missed gene in test (present in ref, not in test)
# ---------------------------------------------------------------------------

@pytest.fixture()
def missed_gene_ref_gff3(tmp_path) -> Path:
    """Reference with two genes."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="gene1", tx_id="tx1",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
        _gene_block(
            chrom="chr1", gene_start=3000, gene_end=4000, strand="-",
            gene_id="gene2", tx_id="tx2",
            exons=((3000, 3400), (3700, 4000)),
            cds=((3050, 3400), (3700, 3950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "ref_missed.gff3", content)


@pytest.fixture()
def missed_gene_test_gff3(tmp_path) -> Path:
    """Test GFF3 with only the first of the two reference genes."""
    content = "\n".join([
        _header(),
        _gene_block(
            gene_id="ENSG00000000001", tx_id="ENST00000000001",
            exons=((1000, 1200), (1600, 2000)),
            cds=((1050, 1200), (1600, 1950)),
        ),
    ]) + "\n"
    return _write_gff3(tmp_path / "test_missed.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: pseudogene recall
# ---------------------------------------------------------------------------

@pytest.fixture()
def pseudo_ref_gff3(tmp_path) -> Path:
    """Reference with one protein-coding gene and one pseudogene."""
    pc = _gene_block(
        gene_id="gene1", tx_id="tx1",
        biotype="protein_coding",
        exons=((1000, 1200), (1600, 2000)),
        cds=((1050, 1200), (1600, 1950)),
    )
    pseudo = _gene_block(
        chrom="chr1", gene_start=5000, gene_end=5500, strand="+",
        gene_id="pseudo1", tx_id="ptx1",
        biotype="pseudogene",
        exons=((5000, 5500),), cds=[(5000, 5500)],
        is_pseudogene=True,
    )
    content = "\n".join([_header(), pc, pseudo]) + "\n"
    return _write_gff3(tmp_path / "ref_pseudo.gff3", content)


@pytest.fixture()
def pseudo_missed_test_gff3(tmp_path) -> Path:
    """Test GFF3 where pseudogene is NOT flagged (regression: called protein_coding)."""
    pc = _gene_block(
        gene_id="ENSG00000000001", tx_id="ENST00000000001",
        biotype="protein_coding",
        exons=((1000, 1200), (1600, 2000)),
        cds=((1050, 1200), (1600, 1950)),
    )
    # Same locus but biotype=protein_coding instead of pseudogene
    pc2 = _gene_block(
        chrom="chr1", gene_start=5000, gene_end=5500, strand="+",
        gene_id="ENSG00000000002", tx_id="ENST00000000002",
        biotype="protein_coding",
        exons=((5000, 5500),), cds=[(5000, 5500)],
        is_pseudogene=False,
    )
    content = "\n".join([_header(), pc, pc2]) + "\n"
    return _write_gff3(tmp_path / "test_pseudo_missed.gff3", content)


# ---------------------------------------------------------------------------
# Fixtures: stats-only GFF3 (multi-biotype, multi-transcript)
# ---------------------------------------------------------------------------

@pytest.fixture()
def stats_gff3(tmp_path) -> Path:
    """
    GFF3 with mixed biotypes for testing compute_stats():
      - 2 protein_coding genes (one has UTR, one doesn't)
      - 1 lnc_RNA gene
      - 1 pseudogene gene
    Gene1 has 2 transcripts (canonical + alternative).
    """
    gene1_pc = "\n".join([
        "##gff-version 3",
        # Gene 1 (protein_coding, 2 transcripts, canonical has UTR)
        "chr1\t.\tgene\t1000\t2500\t.\t+\t.\tID=g1;biotype=protein_coding",
        "chr1\t.\tmRNA\t1000\t2500\t.\t+\t.\tID=t1;Parent=g1;biotype=protein_coding;canonical_transcript=1",
        "chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=t1.e1;Parent=t1",
        "chr1\t.\texon\t1600\t2500\t.\t+\t.\tID=t1.e2;Parent=t1",
        "chr1\t.\tCDS\t1050\t1200\t.\t+\t0\tID=t1.c1;Parent=t1",
        "chr1\t.\tCDS\t1600\t2400\t.\t+\t0\tID=t1.c2;Parent=t1",
        "chr1\t.\tfive_prime_UTR\t1000\t1049\t.\t+\t.\tID=t1.u5;Parent=t1",
        "chr1\t.\tthree_prime_UTR\t2401\t2500\t.\t+\t.\tID=t1.u3;Parent=t1",
        # Alt transcript (no UTR)
        "chr1\t.\tmRNA\t1050\t2400\t.\t+\t.\tID=t1b;Parent=g1;biotype=protein_coding",
        "chr1\t.\texon\t1050\t1200\t.\t+\t.\tID=t1b.e1;Parent=t1b",
        "chr1\t.\texon\t1600\t2400\t.\t+\t.\tID=t1b.e2;Parent=t1b",
        "chr1\t.\tCDS\t1050\t1200\t.\t+\t0\tID=t1b.c1;Parent=t1b",
        "chr1\t.\tCDS\t1600\t2400\t.\t+\t0\tID=t1b.c2;Parent=t1b",
        # Gene 2 (protein_coding, 1 transcript, no UTR)
        "chr1\t.\tgene\t4000\t5000\t.\t-\t.\tID=g2;biotype=protein_coding",
        "chr1\t.\tmRNA\t4000\t5000\t.\t-\t.\tID=t2;Parent=g2;biotype=protein_coding;canonical_transcript=1",
        "chr1\t.\texon\t4000\t4400\t.\t-\t.\tID=t2.e1;Parent=t2",
        "chr1\t.\texon\t4700\t5000\t.\t-\t.\tID=t2.e2;Parent=t2",
        "chr1\t.\tCDS\t4050\t4400\t.\t-\t0\tID=t2.c1;Parent=t2",
        "chr1\t.\tCDS\t4700\t4950\t.\t-\t0\tID=t2.c2;Parent=t2",
        # Gene 3 (lnc_RNA)
        "chr1\t.\tgene\t7000\t8000\t.\t+\t.\tID=g3;biotype=lncRNA",
        "chr1\t.\ttranscript\t7000\t8000\t.\t+\t.\tID=t3;Parent=g3;biotype=lncRNA;canonical_transcript=1",
        "chr1\t.\texon\t7000\t7500\t.\t+\t.\tID=t3.e1;Parent=t3",
        "chr1\t.\texon\t7700\t8000\t.\t+\t.\tID=t3.e2;Parent=t3",
        # Gene 4 (pseudogene)
        "chr1\t.\tgene\t10000\t10500\t.\t+\t.\tID=g4;biotype=pseudogene",
        "chr1\t.\tpseudogenic_transcript\t10000\t10500\t.\t+\t.\tID=t4;Parent=g4;biotype=pseudogene;canonical_transcript=1",
        "chr1\t.\texon\t10000\t10500\t.\t+\t.\tID=t4.e1;Parent=t4",
    ])
    path = tmp_path / "stats_test.gff3"
    path.write_text(gene1_pc + "\n")
    return path
