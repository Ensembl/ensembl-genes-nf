"""
Tests for parse_refseq_gff3.py
"""

import gzip
import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from parse_refseq_gff3 import (
    _infer_biotype,
    _parse_attrs,
    parse_gff3,
    load_synonyms,
    write_gff3,
    GffRecord,
)


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

NCBI_GFF3_BASIC = textwrap.dedent("""\
    ##gff-version 3
    #!gff-spec-version 1.21
    NC_000001.11\tRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene-DDX11L1;Name=DDX11L1;gene_biotype=lncRNA
    NC_000001.11\tRefSeq\tlnc_RNA\t11874\t14409\t.\t+\t.\tID=rna-NR_046018.2;Parent=gene-DDX11L1;Name=NR_046018.2;transcript_biotype=lncRNA
    NC_000001.11\tRefSeq\texon\t11874\t12227\t.\t+\t.\tParent=rna-NR_046018.2
    NC_000001.11\tRefSeq\texon\t12613\t12721\t.\t+\t.\tParent=rna-NR_046018.2
    NC_000001.11\tRefSeq\texon\t13221\t14409\t.\t+\t.\tParent=rna-NR_046018.2
    NC_000001.11\tRefSeq\tgene\t65419\t71585\t.\t+\t.\tID=gene-OR4F5;Name=OR4F5;gene_biotype=protein_coding
    NC_000001.11\tRefSeq\tmRNA\t65419\t71585\t.\t+\t.\tID=rna-NM_001005484.2;Parent=gene-OR4F5;Name=NM_001005484.2
    NC_000001.11\tRefSeq\texon\t65419\t65433\t.\t+\t.\tParent=rna-NM_001005484.2
    NC_000001.11\tRefSeq\texon\t69037\t71585\t.\t+\t.\tParent=rna-NM_001005484.2
    NC_000001.11\tRefSeq\tCDS\t65565\t65433\t.\t+\t0\tParent=rna-NM_001005484.2
    NW_009646201.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-PATCH1;Name=PATCH1;gene_biotype=protein_coding
""")

NCBI_GFF3_NCRNA = textwrap.dedent("""\
    ##gff-version 3
    NC_000001.11\tRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-MIR21;Name=MIR21;gene_biotype=miRNA
    NC_000001.11\tRefSeq\tmiRNA\t1000\t2000\t.\t+\t.\tID=rna-MI0000077;Parent=gene-MIR21;Name=MI0000077
    NC_000001.11\tRefSeq\texon\t1000\t2000\t.\t+\t.\tParent=rna-MI0000077
    NC_000001.11\tRefSeq\tgene\t3000\t4000\t.\t-\t.\tID=gene-SNORD1;Name=SNORD1;gene_biotype=snoRNA
    NC_000001.11\tRefSeq\tsnoRNA\t3000\t4000\t.\t-\t.\tID=rna-SNORD1.1;Parent=gene-SNORD1;Name=SNORD1
    NC_000001.11\tRefSeq\texon\t3000\t4000\t.\t-\t.\tParent=rna-SNORD1.1
""")

SYNONYMS_TSV = "NC_000001.11\t1\nNC_000002.12\t2\n"


def _write_tmp(content: str, suffix: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _write_gz(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(suffix='.gff.gz', delete=False)
    tmp.close()
    with gzip.open(tmp.name, 'wt') as fh:
        fh.write(content)
    return tmp.name


# ---------------------------------------------------------------------------
# _infer_biotype tests
# ---------------------------------------------------------------------------

class TestInferBiotype:
    def test_mrna_is_protein_coding(self):
        assert _infer_biotype('mRNA', {}) == 'protein_coding'

    def test_gene_biotype_overrides_feature(self):
        assert _infer_biotype('gene', {'gene_biotype': 'lncRNA'}) == 'lncRNA'

    def test_lncrna_feature(self):
        assert _infer_biotype('lnc_RNA', {}) == 'lncRNA'

    def test_mirna(self):
        assert _infer_biotype('miRNA', {}) == 'miRNA'

    def test_snorna(self):
        assert _infer_biotype('snoRNA', {}) == 'snoRNA'

    def test_snrna(self):
        assert _infer_biotype('snRNA', {}) == 'snRNA'

    def test_rrna(self):
        assert _infer_biotype('rRNA', {}) == 'rRNA'

    def test_trna(self):
        assert _infer_biotype('tRNA', {}) == 'tRNA'

    def test_unknown_feature_defaults_misc_rna(self):
        assert _infer_biotype('some_novel_feature', {}) == 'misc_RNA'

    def test_transcript_biotype_preferred(self):
        attrs = {'gene_biotype': 'protein_coding', 'transcript_biotype': 'lncRNA'}
        assert _infer_biotype('mRNA', attrs) == 'lncRNA'

    def test_pseudogene(self):
        assert _infer_biotype('gene', {'gene_biotype': 'pseudogene'}) == 'pseudogene'

    def test_lnc_rna_biotype_str(self):
        assert _infer_biotype('gene', {'gene_biotype': 'lnc_RNA'}) == 'lncRNA'


# ---------------------------------------------------------------------------
# _parse_attrs tests
# ---------------------------------------------------------------------------

class TestParseAttrs:
    def test_basic(self):
        attrs = _parse_attrs('ID=gene-DDX11L1;Name=DDX11L1;gene_biotype=lncRNA')
        assert attrs['ID'] == 'gene-DDX11L1'
        assert attrs['Name'] == 'DDX11L1'
        assert attrs['gene_biotype'] == 'lncRNA'

    def test_empty(self):
        assert _parse_attrs('') == {}

    def test_single(self):
        assert _parse_attrs('ID=foo') == {'ID': 'foo'}

    def test_url_encoded_values_preserved(self):
        attrs = _parse_attrs('product=hypothetical%20protein')
        assert attrs['product'] == 'hypothetical%20protein'


# ---------------------------------------------------------------------------
# parse_gff3 tests
# ---------------------------------------------------------------------------

class TestParseGff3:
    def test_parses_genes(self):
        path = _write_tmp(NCBI_GFF3_BASIC, '.gff')
        genes, txs, exons = parse_gff3(path)
        assert len(genes) == 3  # DDX11L1 + OR4F5 + PATCH1
        os.unlink(path)

    def test_parses_transcripts(self):
        path = _write_tmp(NCBI_GFF3_BASIC, '.gff')
        genes, txs, exons = parse_gff3(path)
        assert len(txs) == 2  # NR_046018.2 + NM_001005484.2
        os.unlink(path)

    def test_exons_attached_to_transcript(self):
        path = _write_tmp(NCBI_GFF3_BASIC, '.gff')
        genes, txs, exons = parse_gff3(path)
        assert len(exons['rna-NR_046018.2']) == 3
        os.unlink(path)

    def test_cds_not_counted_as_gene(self):
        path = _write_tmp(NCBI_GFF3_BASIC, '.gff')
        genes, txs, exons = parse_gff3(path)
        # CDS feature should be in exons (treated as exon-like)
        # The gene count should still be 3
        assert len(genes) == 3
        os.unlink(path)

    def test_gz_support(self):
        path = _write_gz(NCBI_GFF3_BASIC)
        genes, txs, exons = parse_gff3(path)
        assert len(genes) == 3
        os.unlink(path)

    def test_gene_record_attributes(self):
        path = _write_tmp(NCBI_GFF3_BASIC, '.gff')
        genes, txs, exons = parse_gff3(path)
        gene = genes.get('gene-DDX11L1')
        assert gene is not None
        assert gene.seqname == 'NC_000001.11'
        assert gene.start == 11874
        assert gene.end == 14409
        assert gene.strand == '+'
        assert gene.attrs['gene_biotype'] == 'lncRNA'
        os.unlink(path)

    def test_ncRNA_types_parsed(self):
        path = _write_tmp(NCBI_GFF3_NCRNA, '.gff')
        genes, txs, exons = parse_gff3(path)
        assert len(genes) == 2
        assert len(txs) == 2
        os.unlink(path)


# ---------------------------------------------------------------------------
# load_synonyms tests
# ---------------------------------------------------------------------------

class TestLoadSynonyms:
    def test_basic_mapping(self):
        path = _write_tmp(SYNONYMS_TSV, '.tsv')
        syn = load_synonyms(path)
        assert syn['NC_000001.11'] == '1'
        assert syn['NC_000002.12'] == '2'
        os.unlink(path)

    def test_none_returns_empty(self):
        assert load_synonyms(None) == {}

    def test_comment_lines_skipped(self):
        path = _write_tmp('# comment\nNC_000001.11\t1\n', '.tsv')
        syn = load_synonyms(path)
        assert 'NC_000001.11' in syn
        os.unlink(path)


# ---------------------------------------------------------------------------
# write_gff3 tests
# ---------------------------------------------------------------------------

class TestWriteGff3:
    def _run(self, gff3_content, synonyms=None, keep_patches=False):
        path = _write_tmp(gff3_content, '.gff')
        genes, txs, exons = parse_gff3(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.gff3', delete=False)
        out.close()
        syn = synonyms or {}
        gene_count, tx_count = write_gff3(genes, txs, exons, out.name, syn, keep_patches)
        with open(out.name) as fh:
            content = fh.read()
        os.unlink(out.name)
        return content, gene_count, tx_count

    def test_writes_header(self):
        content, _, _ = self._run(NCBI_GFF3_BASIC)
        assert '##gff-version 3' in content

    def test_gene_count(self):
        # default: patches excluded → 2 genes (DDX11L1 + OR4F5)
        _, gene_count, _ = self._run(NCBI_GFF3_BASIC)
        assert gene_count == 2

    def test_patch_excluded_by_default(self):
        content, gene_count, _ = self._run(NCBI_GFF3_BASIC)
        assert 'PATCH1' not in content

    def test_keep_patches_includes_nw(self):
        content, gene_count, _ = self._run(NCBI_GFF3_BASIC, keep_patches=True)
        assert gene_count == 3

    def test_synonym_mapping_applied(self):
        syn = {'NC_000001.11': '1'}
        content, _, _ = self._run(NCBI_GFF3_BASIC, synonyms=syn)
        lines = [l for l in content.split('\n') if '\t' in l]
        seqnames = {l.split('\t')[0] for l in lines}
        assert '1' in seqnames
        assert 'NC_000001.11' not in seqnames

    def test_biotype_set_on_gene(self):
        content, _, _ = self._run(NCBI_GFF3_BASIC)
        assert 'biotype=lncRNA' in content
        assert 'biotype=protein_coding' in content

    def test_biotype_best_targeted_replaced_to_refseq(self):
        content, _, _ = self._run(NCBI_GFF3_NCRNA)
        assert 'biotype=miRNA' in content
        assert 'biotype=snoRNA' in content

    def test_exon_features_present(self):
        content, _, _ = self._run(NCBI_GFF3_BASIC)
        lines = [l for l in content.split('\n') if '\t' in l]
        features = {l.split('\t')[2] for l in lines}
        assert 'exon' in features

    def test_transcript_parent_links_gene(self):
        content, _, _ = self._run(NCBI_GFF3_BASIC)
        tx_lines = [l for l in content.split('\n') if '\ttranscript\t' in l]
        for line in tx_lines:
            assert 'Parent=' in line

    def test_transcript_count(self):
        _, gene_count, tx_count = self._run(NCBI_GFF3_BASIC)
        assert tx_count == 2
