#!/usr/bin/env python3
"""
select_canonical.py — Select the canonical transcript for each gene.

Selection criteria (in order):
  1. Longest CDS (sum of CDS feature lengths in bp)
  2. Tie-break: longest total transcript span (end - start + 1)
  3. Tie-break: first in file order

All transcripts are retained in the output.  The canonical transcript gains
the attribute  canonical_transcript=1; all others receive canonical_transcript=0.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


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
# Canonical selection
# ---------------------------------------------------------------------------

def cds_length(tx_id: str, children: Dict[str, List[GffRecord]]) -> int:
    return sum(
        r.end - r.start + 1
        for r in children.get(tx_id, [])
        if r.feature.upper() == 'CDS'
    )


def transcript_span(tx_rec: GffRecord) -> int:
    return tx_rec.end - tx_rec.start + 1


def select_canonical_tx(
    tx_ids: List[str],
    transcripts: Dict[str, GffRecord],
    children: Dict[str, List[GffRecord]],
) -> Optional[str]:
    """Return the ID of the canonical transcript (None if no transcripts)."""
    if not tx_ids:
        return None

    def sort_key(txid: str) -> Tuple[int, int, int]:
        cds_bp = cds_length(txid, children)
        span   = transcript_span(transcripts[txid])
        # Negate so that higher values sort first
        # Position in list gives file-order tie-break (stable sort)
        return (-cds_bp, -span)

    ranked = sorted(tx_ids, key=sort_key)
    return ranked[0]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description='Select the canonical transcript for each gene.'
    )
    ap.add_argument('--gff3', required=True)
    ap.add_argument('--out',  required=True)
    args = ap.parse_args()

    headers, genes, transcripts, children, gene_order, tx_order = parse_gff3(args.gff3)

    n_canonical = 0

    for gid in gene_order:
        tx_ids = tx_order.get(gid, [])
        canonical = select_canonical_tx(tx_ids, transcripts, children)
        if canonical is None:
            continue
        n_canonical += 1
        for txid in tx_ids:
            flag = '1' if txid == canonical else '0'
            transcripts[txid].attrs = _set_attr(
                transcripts[txid].attrs, 'canonical_transcript', flag
            )

    with open(args.out, 'w') as fh:
        for h in headers:
            fh.write(h)
        for gid in gene_order:
            fh.write(_rebuild_raw(genes[gid]) + '\n')
            for txid in tx_order.get(gid, []):
                fh.write(_rebuild_raw(transcripts[txid]) + '\n')
                for child in children.get(txid, []):
                    fh.write(_rebuild_raw(child) + '\n')

    print(f'select_canonical: {n_canonical} canonical transcripts selected', file=sys.stderr)


if __name__ == '__main__':
    main()
