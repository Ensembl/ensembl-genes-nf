"""Round-trip tests for hpp.serialize.

Primary assertion is *encode-stability*: ``encode(decode(encode(x)))`` must equal
``encode(x)``. This proves the serialized form is lossless over everything that
is actually persisted, while remaining robust to deliberately-excluded derived
fields (``source_blocks``, cs-tag caches, spatial indexes). For the core genomic
models we additionally assert full object equality after a round trip.
"""

import pytest

from hpp.models import (
    CDS,
    Exon,
    Gene,
    GenomicInterval,
    Strand,
    SyntenicBlock,
    SyntenicMap,
    Transcript,
)
from hpp.stages.mapping import (
    FeatureMappingResult,
    GeneMappingResult,
    TranscriptMappingResult,
)
from hpp.stages.sex_chromosome import SexChromosomeMap
from hpp.validation.structural import (
    CodonValidation,
    GeneValidation,
    SpliceSiteValidation,
    TranscriptValidation,
)
from hpp.validation.protein import GeneProteinQC, TranscriptProteinQC
from hpp.serialize import decode, encode


def _stable(obj):
    """encode(decode(encode(x))) == encode(x)."""
    e1 = encode(obj)
    e2 = encode(decode(e1))
    return e1 == e2


# --- fixtures -------------------------------------------------------------

def make_interval():
    return GenomicInterval("chr20", 1000, 2000, Strand.MINUS)


def make_exon():
    return Exon(
        feature_id="exon:1",
        feature_type="exon",
        interval=make_interval(),
        attributes={"gene_id": "G1", "tags": ["basic", "ccds"]},
        mapping_identity=0.997,
        exon_number=1,
    )


def make_transcript():
    return Transcript(
        feature_id="tx1",
        feature_type="mRNA",
        interval=make_interval(),
        gene_id="G1",
        biotype="protein_coding",
        attributes={"transcript_id": "tx1"},
        exons=[make_exon()],
        cds_list=[CDS(feature_id="cds1", feature_type="CDS", interval=make_interval(), phase=0)],
    )


def make_gene():
    return Gene(
        feature_id="G1",
        feature_type="gene",
        interval=make_interval(),
        gene_name="GENE1",
        biotype="protein_coding",
        attributes={"gene_id": "G1", "provenance": "projected"},
        transcripts=[make_transcript()],
    )


def make_block():
    return SyntenicBlock(
        ref_interval=make_interval(),
        target_interval=GenomicInterval("20", 1010, 2010, Strand.MINUS),
        identity=0.998,
        alignment_length=1001,
        matches=999,
        mismatches=2,
        cs_tag=":100*ac:50",
        block_id="b1",
    )


# --- core models: full object equality ------------------------------------

@pytest.mark.parametrize("obj", [make_interval(), make_exon(), make_transcript(), make_gene(), make_block()])
def test_core_models_equal_after_roundtrip(obj):
    assert decode(encode(obj)) == obj
    assert _stable(obj)


def test_syntenic_map_blocks_preserved():
    smap = SyntenicMap(blocks=[make_block(), make_block()])
    restored = decode(encode(smap))
    assert restored.blocks == smap.blocks
    assert _stable(smap)


# --- enum, set, tuple primitives ------------------------------------------

def test_strand_enum():
    for s in (Strand.PLUS, Strand.MINUS, Strand.UNSTRANDED):
        assert decode(encode(s)) is s


def test_sex_chromosome_map():
    sx = SexChromosomeMap(
        x_sequences={"chrX", "X"},
        y_sequences={"chrY"},
        par_regions=[("chrX", 10001, 2781479)],
        detection_method="named",
        has_x=True,
        has_y=True,
    )
    restored = decode(encode(sx))
    assert restored == sx
    assert _stable(sx)


# --- mapping result wrappers (source_blocks excluded) ---------------------

def test_mapping_result_stable():
    fmr = FeatureMappingResult(original=make_exon(), mapped=make_exon(), status="mapped")
    tmr = TranscriptMappingResult(
        original=make_transcript(), mapped=make_transcript(), status="mapped",
        exon_results=[fmr], all_exons_mapped=True,
    )
    gmr = GeneMappingResult(
        original=make_gene(), mapped=make_gene(), status="mapped",
        transcript_results=[tmr], transcripts_mapped=1,
    )
    assert _stable(gmr)
    # source_blocks is intentionally dropped -> empty after restore
    restored = decode(encode(fmr))
    assert restored.source_blocks == []


# --- validation + protein qc ----------------------------------------------

def test_validation_results_stable():
    tv = TranscriptValidation(
        transcript_id="tx1",
        splice_sites_checked=3,
        splice_sites_valid=3,
        splice_site_results=[SpliceSiteValidation(0, 100, 200, "GT", "AG", True)],
        start_codon_result=CodonValidation("start", 1000, "ATG", True),
        is_valid=True,
    )
    gv = GeneValidation(gene_id="G1", transcripts_checked=1, transcripts_valid=1,
                        transcript_results=[tv], is_valid=True)
    assert decode(encode(gv)) == gv
    assert _stable(gv)


def test_protein_qc_stable():
    tq = TranscriptProteinQC(transcript_id="tx1", reference_transcript_id="rtx1",
                             status="ok", identity=0.99, coverage=0.98)
    gq = GeneProteinQC(gene_id="G1", mapped_from="rG1", status="ok",
                       transcripts_checked=1, transcripts_ok=1, transcript_results=[tq])
    assert decode(encode(gq)) == gq
    assert _stable(gq)
