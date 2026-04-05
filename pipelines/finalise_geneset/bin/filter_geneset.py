#!/usr/bin/env python3
"""
filter_geneset.py — Remove poor-quality transcripts from a gene GFF3.

Filters applied:
  1. Short ORF: transcripts whose total CDS length encodes fewer than
     --min-orf-aa amino acids are removed.
  2. Tiny intron: transcripts with any intron shorter than --min-intron-size bp
     are removed (frameshifted model indicator).
  3. If all transcripts of a gene are removed, the gene is also removed.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GffRecord:
    """One line of a GFF3 file (excluding comment/directive lines)."""
    seqname: str
    source:  str
    feature: str
    start:   int
    end:     int
    score:   str
    strand:  str
    phase:   str
    attrs:   str
    raw:     str  # original line (without newline) — used for re-emission


def _parse_attrs(attrs: str) -> Dict[str, str]:
    """Return key→value dict from a GFF3 attribute column."""
    d: Dict[str, str] = {}
    for token in attrs.rstrip(';').split(';'):
        token = token.strip()
        if '=' in token:
            k, v = token.split('=', 1)
            d[k] = v
    return d


def _set_attr(attrs: str, key: str, value: str) -> str:
    """Set (or add) a single key=value pair in a GFF3 attribute string."""
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


def parse_gff3(path: str):
    """
    Parse a GFF3 file and return:
      - headers: list of comment/directive lines (strings, with newline)
      - genes:   dict gene_id → GffRecord (gene feature)
      - transcripts: dict tx_id → GffRecord (mRNA/transcript feature)
      - children: dict tx_id → list of GffRecord (exon, CDS, UTR, …)
      - gene_order: list of gene_ids in file order
      - tx_order: dict gene_id → list of tx_ids in file order
    """
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
                seqname=cols[0],
                source=cols[1],
                feature=cols[2],
                start=int(cols[3]),
                end=int(cols[4]),
                score=cols[5],
                strand=cols[6],
                phase=cols[7],
                attrs=cols[8],
                raw=raw,
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
                # exon, CDS, UTR, etc.
                for parent in a.get('Parent', '').split(','):
                    parent = parent.strip()
                    if parent:
                        children[parent].append(rec)

    return headers, genes, transcripts, children, gene_order, tx_order


# ---------------------------------------------------------------------------
# Filtering logic
# ---------------------------------------------------------------------------

def cds_length(tx_id: str, children: Dict[str, List[GffRecord]]) -> int:
    """Total CDS length (bp) for a transcript."""
    total = 0
    for rec in children.get(tx_id, []):
        if rec.feature.upper() == 'CDS':
            total += rec.end - rec.start + 1
    return total


def min_intron_size_for_tx(tx_id: str, children: Dict[str, List[GffRecord]]) -> Optional[int]:
    """
    Return the size of the smallest intron implied by the exon structure of
    this transcript, or None if there is only one exon (no introns).
    """
    exons = sorted(
        [r for r in children.get(tx_id, []) if r.feature.lower() == 'exon'],
        key=lambda r: r.start,
    )
    if len(exons) < 2:
        return None
    intron_sizes = []
    for i in range(1, len(exons)):
        intron_size = exons[i].start - exons[i - 1].end - 1
        intron_sizes.append(intron_size)
    return min(intron_sizes) if intron_sizes else None


def should_remove_transcript(
    tx_id: str,
    children: Dict[str, List[GffRecord]],
    min_orf_aa: int,
    min_intron_size: int,
) -> Tuple[bool, str]:
    """
    Return (True, reason) if the transcript should be removed, else (False, '').
    """
    cds_bp = cds_length(tx_id, children)
    if cds_bp < min_orf_aa * 3:
        return True, f'short_orf cds_bp={cds_bp} < {min_orf_aa * 3}'

    smallest = min_intron_size_for_tx(tx_id, children)
    if smallest is not None and smallest < min_intron_size:
        return True, f'tiny_intron size={smallest} < {min_intron_size}'

    return False, ''


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def _rebuild_raw(rec: GffRecord) -> str:
    """Reconstruct the tab-separated GFF3 line from a (possibly modified) record."""
    return '\t'.join([
        rec.seqname, rec.source, rec.feature,
        str(rec.start), str(rec.end),
        rec.score, rec.strand, rec.phase,
        rec.attrs,
    ])


def write_filtered_gff3(
    out_path: str,
    headers: List[str],
    genes: Dict[str, GffRecord],
    transcripts: Dict[str, GffRecord],
    children: Dict[str, List[GffRecord]],
    gene_order: List[str],
    tx_order: Dict[str, List[str]],
    removed_txs: set,
    removed_genes: set,
) -> None:
    with open(out_path, 'w') as fh:
        for h in headers:
            fh.write(h)
        for gid in gene_order:
            if gid in removed_genes:
                continue
            gene_rec = genes[gid]
            fh.write(_rebuild_raw(gene_rec) + '\n')
            for txid in tx_order.get(gid, []):
                if txid in removed_txs:
                    continue
                tx_rec = transcripts[txid]
                fh.write(_rebuild_raw(tx_rec) + '\n')
                for child in children.get(txid, []):
                    fh.write(_rebuild_raw(child) + '\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description='Filter poor-quality transcripts from a gene GFF3.'
    )
    ap.add_argument('--gff3',           required=True, help='Input GFF3 file')
    ap.add_argument('--out',            required=True, help='Output GFF3 file')
    ap.add_argument('--min-orf-aa',     type=int, default=100,
                    help='Minimum CDS length in amino acids (default: 100)')
    ap.add_argument('--min-intron-size', type=int, default=10,
                    help='Minimum intron size in bp (default: 10)')
    args = ap.parse_args()

    headers, genes, transcripts, children, gene_order, tx_order = parse_gff3(args.gff3)

    removed_txs:   set = set()
    removed_genes: set = set()

    n_tx_removed = 0
    n_gene_removed = 0
    reasons: Dict[str, int] = defaultdict(int)

    for gid in gene_order:
        kept_txs = []
        for txid in tx_order.get(gid, []):
            remove, reason = should_remove_transcript(
                txid, children, args.min_orf_aa, args.min_intron_size
            )
            if remove:
                removed_txs.add(txid)
                n_tx_removed += 1
                reasons[reason.split()[0]] += 1
            else:
                kept_txs.append(txid)

        if not kept_txs:
            removed_genes.add(gid)
            n_gene_removed += 1

    write_filtered_gff3(
        args.out,
        headers,
        genes,
        transcripts,
        children,
        gene_order,
        tx_order,
        removed_txs,
        removed_genes,
    )

    total_tx = sum(len(v) for v in tx_order.values())
    total_genes = len(genes)
    print(
        f'filter_geneset: {total_genes} genes, {total_tx} transcripts in; '
        f'{n_tx_removed} transcripts removed ({reasons}); '
        f'{n_gene_removed} genes removed',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
