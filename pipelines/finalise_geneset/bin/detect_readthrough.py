#!/usr/bin/env python3
"""
detect_readthrough.py — Flag transcripts that span two independent gene loci.

A readthrough transcript is one whose exon coordinates overlap with CDS exons
of TWO OR MORE different genes (other than its own parent gene) on the same
strand of the same chromosome.

Flagging:
  - transcript biotype → 'readthrough_transcript'
  - gene biotype → 'readthrough' if ALL its transcripts are flagged
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class GffRecord:
    seqname: str
    source:  str
    feature: str
    start:   int
    end:     int
    score:   str
    strand:  str
    phase:   str
    attrs:   str
    raw:     str


def _parse_attrs(attrs: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for token in attrs.rstrip(';').split(';'):
        token = token.strip()
        if '=' in token:
            k, v = token.split('=', 1)
            d[k] = v
    return d


def _set_attr(attrs: str, key: str, value: str) -> str:
    pairs = []
    found = False
    for token in attrs.rstrip(';').split(';'):
        token = token.strip()
        if not token:
            continue
        if '=' in token:
            k, _ = token.split('=', 1)
            if k == key:
                pairs.append(f'{key}={value}')
                found = True
                continue
        pairs.append(token)
    if not found:
        pairs.append(f'{key}={value}')
    return ';'.join(pairs)


def _rebuild_raw(rec: GffRecord) -> str:
    return '\t'.join([
        rec.seqname, rec.source, rec.feature,
        str(rec.start), str(rec.end),
        rec.score, rec.strand, rec.phase,
        rec.attrs,
    ])


def parse_gff3(path: str):
    headers = []
    genes: Dict[str, GffRecord] = {}
    transcripts: Dict[str, GffRecord] = {}
    children: Dict[str, List[GffRecord]] = defaultdict(list)
    gene_order: List[str] = []
    tx_order: Dict[str, List[str]] = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            raw = line.rstrip('\n')
            if raw.startswith('#') or raw.strip() == '':
                headers.append(raw + '\n')
                continue
            cols = raw.split('\t')
            if len(cols) < 9:
                continue
            rec = GffRecord(
                seqname=cols[0], source=cols[1], feature=cols[2],
                start=int(cols[3]), end=int(cols[4]),
                score=cols[5], strand=cols[6], phase=cols[7],
                attrs=cols[8], raw=raw,
            )
            a = _parse_attrs(rec.attrs)
            feat = rec.feature.lower()
            if feat == 'gene':
                gid = a.get('ID', '')
                genes[gid] = rec
                gene_order.append(gid)
            elif feat in ('mrna', 'transcript'):
                txid = a.get('ID', '')
                gid  = a.get('Parent', '')
                transcripts[txid] = rec
                tx_order[gid].append(txid)
            else:
                for parent in a.get('Parent', '').split(','):
                    parent = parent.strip()
                    if parent:
                        children[parent].append(rec)

    return headers, genes, transcripts, children, gene_order, tx_order


# ---------------------------------------------------------------------------
# Build a CDS interval index per (seqname, strand, gene_id)
# ---------------------------------------------------------------------------

def build_cds_index(
    gene_order: List[str],
    genes: Dict[str, GffRecord],
    tx_order: Dict[str, List[str]],
    children: Dict[str, List[GffRecord]],
) -> Dict[Tuple[str, str], List[Tuple[int, int, str]]]:
    """
    Return a dict (seqname, strand) → sorted list of (start, end, gene_id)
    for every CDS feature in the GFF3.
    """
    index: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
    for gid in gene_order:
        gene_rec = genes[gid]
        key = (gene_rec.seqname, gene_rec.strand)
        for txid in tx_order.get(gid, []):
            for child in children.get(txid, []):
                if child.feature.upper() == 'CDS':
                    index[key].append((child.start, child.end, gid))
    for key in index:
        index[key].sort()
    return index


def overlapping_genes(
    exon_start: int,
    exon_end: int,
    cds_list: List[Tuple[int, int, str]],
    own_gene_id: str,
) -> Set[str]:
    """Return set of gene IDs (excluding own_gene_id) whose CDS overlaps [exon_start, exon_end]."""
    result: Set[str] = set()
    for (cs, ce, gid) in cds_list:
        if cs > exon_end:
            break
        if ce < exon_start:
            continue
        if gid != own_gene_id:
            result.add(gid)
    return result


def is_readthrough(
    tx_id: str,
    parent_gene_id: str,
    tx_rec: GffRecord,
    children: Dict[str, List[GffRecord]],
    cds_index: Dict[Tuple[str, str], List[Tuple[int, int, str]]],
) -> bool:
    """
    A transcript is readthrough if its exons overlap the CDS of 2+ genes
    other than its own parent gene.
    """
    key = (tx_rec.seqname, tx_rec.strand)
    cds_list = cds_index.get(key, [])
    if not cds_list:
        return False

    hit_genes: Set[str] = set()
    for child in children.get(tx_id, []):
        if child.feature.lower() == 'exon':
            hit_genes |= overlapping_genes(child.start, child.end, cds_list, parent_gene_id)
        if len(hit_genes) >= 2:
            return True
    return len(hit_genes) >= 2


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description='Flag readthrough transcripts in a gene GFF3.'
    )
    ap.add_argument('--gff3',    required=True)
    ap.add_argument('--out',     required=True)
    ap.add_argument('--max-gap', type=int, default=10000,
                    help='Max intergenic gap to consider for readthrough (default: 10000)')
    args = ap.parse_args()

    headers, genes, transcripts, children, gene_order, tx_order = parse_gff3(args.gff3)
    cds_index = build_cds_index(gene_order, genes, tx_order, children)

    n_rt_tx    = 0
    n_rt_gene  = 0

    for gid in gene_order:
        tx_ids = tx_order.get(gid, [])
        if not tx_ids:
            continue

        rt_flags = {}
        for txid in tx_ids:
            tx_rec = transcripts[txid]
            rt = is_readthrough(txid, gid, tx_rec, children, cds_index)
            rt_flags[txid] = rt
            if rt:
                transcripts[txid].attrs = _set_attr(
                    transcripts[txid].attrs, 'biotype', 'readthrough_transcript'
                )
                n_rt_tx += 1

        if all(rt_flags.values()) and any(rt_flags.values()):
            genes[gid].attrs = _set_attr(genes[gid].attrs, 'biotype', 'readthrough')
            n_rt_gene += 1

    with open(args.out, 'w') as fh:
        for h in headers:
            fh.write(h)
        for gid in gene_order:
            fh.write(_rebuild_raw(genes[gid]) + '\n')
            for txid in tx_order.get(gid, []):
                fh.write(_rebuild_raw(transcripts[txid]) + '\n')
                for child in children.get(txid, []):
                    fh.write(_rebuild_raw(child) + '\n')

    print(
        f'detect_readthrough: {n_rt_tx} readthrough transcripts, '
        f'{n_rt_gene} readthrough genes flagged',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
