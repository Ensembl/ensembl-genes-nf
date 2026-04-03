"""
Tests for parse_assembly_report.py
"""

import gzip
import json
import os
import sys
import tempfile
import textwrap

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from parse_assembly_report import (
    parse_assembly_report,
    write_synonyms,
    write_metadata,
)


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

REPORT_BASIC = textwrap.dedent("""\
    # Assembly name:  GRCh38.p14
    # Organism name:  Homo sapiens (human)
    # Taxid:          9606
    # Assembly accession: GCA_000001405.29
    # RefSeq assembly accession: GCF_000001405.40
    # Assembly level: Chromosome
    # Genome representation: full
    # Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
    1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1
    2\tassembled-molecule\t2\tChromosome\tCM000664.1\t=\tNC_000002.12\tPrimary Assembly\t242193529\tchr2
    X\tassembled-molecule\tX\tChromosome\tCM000685.2\t=\tNC_000023.11\tPrimary Assembly\t156040895\tchrX
    MT\tassembled-molecule\tMT\tMitochondrion\tJ01415.2\t=\tNC_012920.1\tNon-nuclear\t16569\tchrM
    HSCHR1_CTG1_UNLOCALIZED\tunlocalized-scaffold\t1\tChromosome\tGL000006.2\t=\tNT_187361.1\tPrimary Assembly\t171115067\tchr1_GL000006v2_random
    SCAFFOLD123\tunplaced-scaffold\tna\tna\tKV880768.2\t=\tNW_019805495.1\tPrimary Assembly\t182439\tna
""")

REPORT_NO_REFSEQ = textwrap.dedent("""\
    # Assembly name:  CustomAssembly
    # Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
    chr1\tassembled-molecule\t1\tChromosome\tCM000001.1\t=\tna\tPrimary Assembly\t1000000\tna
""")


def _write_tmp(content: str, suffix='.txt') -> str:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _write_gz(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(suffix='.txt.gz', delete=False)
    tmp.close()
    with gzip.open(tmp.name, 'wt') as fh:
        fh.write(content)
    return tmp.name


# ---------------------------------------------------------------------------
# parse_assembly_report tests
# ---------------------------------------------------------------------------

class TestParseAssemblyReport:
    def test_parses_all_sequences(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        assert len(rows) == 6
        os.unlink(path)

    def test_metadata_extracted(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        assert meta.get('assembly_name') == 'GRCh38.p14'
        assert meta.get('taxid') == '9606'
        os.unlink(path)

    def test_sequence_row_fields(self):
        path = _write_tmp(REPORT_BASIC)
        rows, _ = parse_assembly_report(path)
        chr1 = next(r for r in rows if r['seq_name'] == '1')
        assert chr1['refseq'] == 'NC_000001.11'
        assert chr1['ucsc'] == 'chr1'
        assert chr1['genbank'] == 'CM000663.2'
        assert chr1['length'] == 248956422
        assert chr1['role'] == 'assembled-molecule'
        os.unlink(path)

    def test_na_refseq_is_none(self):
        path = _write_tmp(REPORT_NO_REFSEQ)
        rows, _ = parse_assembly_report(path)
        assert rows[0]['refseq'] is None
        os.unlink(path)

    def test_na_ucsc_is_none(self):
        path = _write_tmp(REPORT_BASIC)
        rows, _ = parse_assembly_report(path)
        scaffold = next(r for r in rows if r['seq_name'] == 'SCAFFOLD123')
        assert scaffold['ucsc'] is None
        os.unlink(path)

    def test_gz_support(self):
        path = _write_gz(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        assert len(rows) == 6
        os.unlink(path)

    def test_header_comment_lines_skipped(self):
        path = _write_tmp(REPORT_BASIC)
        rows, _ = parse_assembly_report(path)
        # The column-header comment line should not produce a row
        assert not any(r['seq_name'].startswith('Sequence') for r in rows)
        os.unlink(path)

    def test_mt_sequence_present(self):
        path = _write_tmp(REPORT_BASIC)
        rows, _ = parse_assembly_report(path)
        mt = next((r for r in rows if r['seq_name'] == 'MT'), None)
        assert mt is not None
        assert mt['refseq'] == 'NC_012920.1'
        os.unlink(path)


# ---------------------------------------------------------------------------
# write_synonyms tests
# ---------------------------------------------------------------------------

class TestWriteSynonyms:
    def _rows(self):
        path = _write_tmp(REPORT_BASIC)
        rows, _ = parse_assembly_report(path)
        os.unlink(path)
        return rows

    def test_skips_rows_without_refseq(self):
        path = _write_tmp(REPORT_NO_REFSEQ)
        rows, _ = parse_assembly_report(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        out.close()
        n = write_synonyms(rows, out.name)
        assert n == 0
        os.unlink(out.name)

    def test_returns_count_with_refseq(self):
        rows = self._rows()
        out = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        out.close()
        n = write_synonyms(rows, out.name)
        # 6 rows, all have refseq in REPORT_BASIC
        assert n == 6
        os.unlink(out.name)

    def test_ucsc_name_preferred_over_seqname(self):
        rows = self._rows()
        out = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        out.close()
        write_synonyms(rows, out.name)
        with open(out.name) as fh:
            content = fh.read()
        # chr1 should appear as the synonym for NC_000001.11
        assert 'NC_000001.11\tchr1\t' in content
        os.unlink(out.name)

    def test_genbank_fallback_when_no_ucsc(self):
        rows = self._rows()
        out = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        out.close()
        write_synonyms(rows, out.name)
        with open(out.name) as fh:
            content = fh.read()
        # SCAFFOLD123 has ucsc=None, should use genbank KV880768.2
        assert 'NW_019805495.1\tKV880768.2\t' in content
        os.unlink(out.name)

    def test_header_line_written(self):
        rows = self._rows()
        out = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
        out.close()
        write_synonyms(rows, out.name)
        with open(out.name) as fh:
            first_line = fh.readline()
        assert first_line.startswith('#')
        os.unlink(out.name)


# ---------------------------------------------------------------------------
# write_metadata tests
# ---------------------------------------------------------------------------

class TestWriteMetadata:
    def test_writes_valid_json(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
        out.close()
        write_metadata(meta, rows, out.name)
        with open(out.name) as fh:
            doc = json.load(fh)
        assert 'assembly_metadata' in doc
        assert 'sequence_summary' in doc
        os.unlink(out.name)

    def test_total_sequence_count(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
        out.close()
        write_metadata(meta, rows, out.name)
        with open(out.name) as fh:
            doc = json.load(fh)
        assert doc['sequence_summary']['total_sequences'] == 6
        os.unlink(out.name)

    def test_total_length_computed(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
        out.close()
        write_metadata(meta, rows, out.name)
        with open(out.name) as fh:
            doc = json.load(fh)
        assert doc['sequence_summary']['total_length'] > 0
        os.unlink(out.name)

    def test_assembly_name_in_meta(self):
        path = _write_tmp(REPORT_BASIC)
        rows, meta = parse_assembly_report(path)
        os.unlink(path)
        out = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
        out.close()
        write_metadata(meta, rows, out.name)
        with open(out.name) as fh:
            doc = json.load(fh)
        assert doc['assembly_metadata'].get('assembly_name') == 'GRCh38.p14'
        os.unlink(out.name)
