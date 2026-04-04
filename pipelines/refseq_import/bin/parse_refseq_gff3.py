#!/usr/bin/env python3
"""
parse_refseq_gff3.py — Parse NCBI RefSeq GFF3 and emit Ensembl-style GFF3.

NCBI RefSeq GFF3 format notes:
  - Sequence names are RefSeq accessions (NC_000001.11, NW_..., etc.)
  - gene_biotype / transcript_biotype attributes specify RNA class
  - Feature types: gene, mRNA, lnc_RNA, miRNA, snoRNA, snRNA, rRNA, tRNA,
    ncRNA, primary_transcript, exon, CDS, transcript, misc_RNA, ...
  - Dbxref contains GeneID, HGNC, MIM, RefSeq etc.
  - Parent links transcript → gene and exon/CDS → transcript

Output:
  - Standard GFF3 with gene/transcript/exon hierarchy
  - biotype= from gene_biotype or inferred from feature type
  - seq names optionally remapped via --synonyms TSV (refseq_acc → chr_name)
  - Non-genomic sequences (NT_, NW_) optionally filtered out

Usage:
    parse_refseq_gff3.py \\
        --gff3     GCF_000001405.40_GRCh38.p14_genomic.gff.gz \\
        --out      refseq.gff3 \\
        [--synonyms synonyms.tsv] \\
        [--keep_patches]          # include NT_/NW_ scaffolds
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Biotype inference
# ---------------------------------------------------------------------------

# Map NCBI feature types / gene_biotype values → Ensembl biotypes
_FEATURE_BIOTYPE: Dict[str, str] = {
    'mRNA':                'protein_coding',
    'lnc_RNA':             'lncRNA',
    'miRNA':               'miRNA',
    'snoRNA':              'snoRNA',
    'snRNA':               'snRNA',
    'rRNA':                'rRNA',
    'tRNA':                'tRNA',
    'ncRNA':               'misc_RNA',
    'misc_RNA':            'misc_RNA',
    'primary_transcript':  'misc_RNA',
    'transcript':          'misc_RNA',
    'vault_RNA':           'vault_RNA',
    'scRNA':               'scaRNA',
    'Y_RNA':               'Y_RNA',
    'SRP_RNA':             'SRP_RNA',
    'RNase_MRP_RNA':       'RNase_MRP_RNA',
    'antisense_RNA':       'antisense',
    'pseudogene':          'pseudogene',
    'processed_pseudogene':'processed_pseudogene',
}

_GENE_BIOTYPE_MAP: Dict[str, str] = {
    'protein_coding':           'protein_coding',
    'lncRNA':                   'lncRNA',
    'lnc_RNA':                  'lncRNA',
    'miRNA':                    'miRNA',
    'snoRNA':                   'snoRNA',
    'snRNA':                    'snRNA',
    'rRNA':                     'rRNA',
    'tRNA':                     'tRNA',
    'misc_RNA':                 'misc_RNA',
    'pseudogene':               'pseudogene',
    'transcribed_pseudogene':   'transcribed_unprocessed_pseudogene',
    'processed_pseudogene':     'processed_pseudogene',
    'unprocessed_pseudogene':   'unprocessed_pseudogene',
    'other':                    'misc_RNA',
    'vault_RNA':                'vault_RNA',
    'Y_RNA':                    'Y_RNA',
    'SRP_RNA':                  'SRP_RNA',
    'scRNA':                    'scaRNA',
    'antisense_RNA':            'antisense',
}


def _infer_biotype(feature: str, attrs: dict) -> str:
    """Derive Ensembl biotype from GFF3 feature type and attributes."""
    # transcript_biotype is most specific; gene_biotype is a fallback
    for key in ('transcript_biotype', 'gene_biotype', 'biotype'):
        val = attrs.get(key, '')
        if val:
            mapped = _GENE_BIOTYPE_MAP.get(val)
            if mapped:
                return mapped
    # Fall back to feature type
    return _FEATURE_BIOTYPE.get(feature, 'misc_RNA')


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

_ATTR_RE = re.compile(r'([^=;]+)=([^;]*)')


def _parse_attrs(attr_str: str) -> dict:
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr_str)}


def _open_gff(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


@dataclass
class GffRecord:
    seqname: str
    source:  str
    feature: str
    start:   int
    end:     int
    score:   str
    strand:  str
    frame:   str
    attrs:   dict
    raw_id:  str = ''
    parent:  str = ''


_GENE_FEATURES = frozenset(['gene', 'pseudogene'])
_TX_FEATURES   = frozenset([
    'mRNA', 'lnc_RNA', 'miRNA', 'snoRNA', 'snRNA', 'rRNA', 'tRNA',
    'ncRNA', 'misc_RNA', 'primary_transcript', 'transcript',
    'vault_RNA', 'scRNA', 'Y_RNA', 'SRP_RNA', 'RNase_MRP_RNA',
    'antisense_RNA', 'processed_pseudogene',
])
_EXON_FEATURES = frozenset(['exon', 'CDS'])


def parse_gff3(path: str) -> Tuple[
    Dict[str, GffRecord],   # id → gene record
    Dict[str, GffRecord],   # id → transcript record
    Dict[str, List[GffRecord]],  # tx_id → list of exon records
]:
    genes: Dict[str, GffRecord] = {}
    txs:   Dict[str, GffRecord] = {}
    exons: Dict[str, List[GffRecord]] = defaultdict(list)

    with _open_gff(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in (_GENE_FEATURES | _TX_FEATURES | _EXON_FEATURES):
                continue

            attrs = _parse_attrs(cols[8])
            rec = GffRecord(
                seqname = cols[0],
                source  = cols[1],
                feature = feature,
                start   = int(cols[3]),
                end     = int(cols[4]),
                score   = cols[5],
                strand  = cols[6],
                frame   = cols[7],
                attrs   = attrs,
                raw_id  = attrs.get('ID', ''),
                parent  = attrs.get('Parent', ''),
            )

            if feature in _GENE_FEATURES:
                genes[rec.raw_id] = rec
            elif feature in _TX_FEATURES:
                txs[rec.raw_id] = rec
            elif feature in _EXON_FEATURES:
                # attach to parent transcript
                if rec.parent:
                    exons[rec.parent].append(rec)

    return genes, txs, exons


# ---------------------------------------------------------------------------
# Synonym mapping (RefSeq accession → chromosome name)
# ---------------------------------------------------------------------------

def load_synonyms(path: Optional[str]) -> Dict[str, str]:
    """Load TSV: col0=refseq_accession, col1=chr_name (e.g. NC_000001.11→1)"""
    if not path:
        return {}
    mapping: Dict[str, str] = {}
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

# Sequences to exclude by default (patches, unlocalized)
_EXCLUDE_RE = re.compile(r'^N[TW]_')


def write_gff3(
    genes:    Dict[str, GffRecord],
    txs:      Dict[str, GffRecord],
    exons:    Dict[str, List[GffRecord]],
    out_path: str,
    synonyms: Dict[str, str],
    keep_patches: bool,
) -> Tuple[int, int]:
    """Write Ensembl-style GFF3.  Returns (gene_count, transcript_count)."""
    gene_count = 0
    tx_count   = 0

    # Build gene→transcript index
    gene_to_txs: Dict[str, List[str]] = defaultdict(list)
    for tx_id, tx in txs.items():
        # Parent may be a comma-separated list (rare); take first
        parent = tx.parent.split(',')[0] if tx.parent else ''
        gene_to_txs[parent].append(tx_id)

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')

        for gene_id, gene in genes.items():
            seqname = synonyms.get(gene.seqname, gene.seqname)

            # Filter out patch sequences unless requested
            if not keep_patches and _EXCLUDE_RE.match(gene.seqname):
                continue

            biotype  = _infer_biotype(gene.feature, gene.attrs)
            gene_name = gene.attrs.get('Name', gene.attrs.get('gene', gene_id))
            gene_ensid = f'refseq_gene_{gene_id}'

            fh.write(
                f'{seqname}\tRefSeq\tgene\t{gene.start}\t{gene.end}\t'
                f'{gene.score}\t{gene.strand}\t.\t'
                f'ID={gene_ensid};Name={gene_name};biotype={biotype}\n'
            )
            gene_count += 1

            child_tx_ids = gene_to_txs.get(gene_id, [])
            for tx_id in child_tx_ids:
                tx = txs[tx_id]
                tx_biotype = _infer_biotype(tx.feature, tx.attrs)
                tx_name    = tx.attrs.get('Name', tx.attrs.get('transcript_id', tx_id))
                tx_ensid   = f'refseq_tx_{tx_id}'

                fh.write(
                    f'{seqname}\tRefSeq\ttranscript\t{tx.start}\t{tx.end}\t'
                    f'{tx.score}\t{tx.strand}\t.\t'
                    f'ID={tx_ensid};Parent={gene_ensid};Name={tx_name};biotype={tx_biotype}\n'
                )
                tx_count += 1

                tx_exons = sorted(exons.get(tx_id, []), key=lambda e: e.start)
                for i, exon in enumerate(tx_exons, 1):
                    fh.write(
                        f'{seqname}\tRefSeq\texon\t{exon.start}\t{exon.end}\t'
                        f'.\t{exon.strand}\t.\t'
                        f'ID={tx_ensid}_exon_{i};Parent={tx_ensid}\n'
                    )

    return gene_count, tx_count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff3',         required=True, help='NCBI RefSeq GFF3 (.gz ok)')
    ap.add_argument('--out',          required=True, help='Output GFF3 path')
    ap.add_argument('--synonyms',     default=None,
                    help='TSV mapping RefSeq accession → sequence name')
    ap.add_argument('--keep_patches', action='store_true',
                    help='Include NT_/NW_ patch/scaffold sequences')
    args = ap.parse_args()

    synonyms = load_synonyms(args.synonyms)
    genes, txs, exons = parse_gff3(args.gff3)
    gene_count, tx_count = write_gff3(genes, txs, exons, args.out,
                                      synonyms, args.keep_patches)
    print(f'parse_refseq_gff3: {gene_count} genes, {tx_count} transcripts written',
          file=sys.stderr)


if __name__ == '__main__':
    main()
