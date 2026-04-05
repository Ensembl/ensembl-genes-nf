#!/usr/bin/env python3
"""
detect_pseudogenes.py — Flag pseudogenes in a gene GFF3 using repeat overlap.

Two criteria:
  1. Single-exon coding gene: if the fraction of CDS bp covered by repeats
     exceeds --max-repeat-coverage, set biotype='processed_pseudogene'.
  2. Multi-exon coding gene: if ALL introns are frameshifted (length < 60 bp),
     set biotype='pseudogene'.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Re-use parsing helpers (duplicated here so each script is self-contained)
# ---------------------------------------------------------------------------

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


def parse_repeats(path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse a repeat GFF3 and return a dict:
      seqname → sorted list of (start, end) intervals (1-based, closed).
    """
    repeats: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                continue
            try:
                start = int(cols[3])
                end   = int(cols[4])
            except ValueError:
                continue
            repeats[cols[0]].append((start, end))
    for chrom in repeats:
        repeats[chrom].sort()
    return repeats


# ---------------------------------------------------------------------------
# Coverage calculation
# ---------------------------------------------------------------------------

def covered_bases(
    intervals: List[Tuple[int, int]],
    query_start: int,
    query_end: int,
) -> int:
    """
    Count bp of [query_start, query_end] covered by any interval in *intervals*
    (already sorted by start).  All coordinates 1-based, closed.
    """
    covered = 0
    for (s, e) in intervals:
        if s > query_end:
            break
        ov_s = max(s, query_start)
        ov_e = min(e, query_end)
        if ov_s <= ov_e:
            covered += ov_e - ov_s + 1
    return covered


def repeat_cds_coverage(
    tx_id: str,
    children: Dict[str, List[GffRecord]],
    repeat_ivs: List[Tuple[int, int]],
) -> float:
    """Return fraction of total CDS bp covered by repeats."""
    cds_recs = [r for r in children.get(tx_id, []) if r.feature.upper() == 'CDS']
    if not cds_recs:
        return 0.0
    total_cds = sum(r.end - r.start + 1 for r in cds_recs)
    cov = sum(
        covered_bases(repeat_ivs, r.start, r.end) for r in cds_recs
    )
    return cov / total_cds if total_cds > 0 else 0.0


# ---------------------------------------------------------------------------
# Intron frameshift detection
# ---------------------------------------------------------------------------

FRAMESHIFT_INTRON_MAX = 60  # introns shorter than this are suspect


def all_introns_frameshifted(
    tx_id: str,
    children: Dict[str, List[GffRecord]],
) -> bool:
    """
    Return True if there are >=1 introns AND every intron is shorter than
    FRAMESHIFT_INTRON_MAX bp.
    """
    exons = sorted(
        [r for r in children.get(tx_id, []) if r.feature.lower() == 'exon'],
        key=lambda r: r.start,
    )
    if len(exons) < 2:
        return False
    intron_sizes = [
        exons[i].start - exons[i - 1].end - 1
        for i in range(1, len(exons))
    ]
    return all(s < FRAMESHIFT_INTRON_MAX for s in intron_sizes)


def exon_count(tx_id: str, children: Dict[str, List[GffRecord]]) -> int:
    return sum(1 for r in children.get(tx_id, []) if r.feature.lower() == 'exon')


def has_cds(tx_id: str, children: Dict[str, List[GffRecord]]) -> bool:
    return any(r.feature.upper() == 'CDS' for r in children.get(tx_id, []))


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def classify_gene(
    gid: str,
    tx_ids: List[str],
    children: Dict[str, List[GffRecord]],
    repeat_ivs: List[Tuple[int, int]],
    max_repeat_coverage: float,
) -> Optional[str]:
    """
    Return 'processed_pseudogene', 'pseudogene', or None.

    Strategy: classify based on the representative transcript (first one with CDS).
    """
    coding_txs = [t for t in tx_ids if has_cds(t, children)]
    if not coding_txs:
        return None

    # Single-exon coding genes — check repeat coverage of CDS
    # (all coding transcripts must be single-exon)
    all_single = all(exon_count(t, children) <= 1 for t in coding_txs)
    if all_single:
        # Use the transcript with the most CDS coverage for the verdict
        rep_tx = coding_txs[0]
        cov = repeat_cds_coverage(rep_tx, children, repeat_ivs)
        if cov > max_repeat_coverage:
            return 'processed_pseudogene'
        return None

    # Multi-exon coding genes — check if ALL coding transcripts have all frameshifted introns
    all_frameshifted = all(
        all_introns_frameshifted(t, children)
        for t in coding_txs
        if exon_count(t, children) > 1
    )
    multi_exon_txs = [t for t in coding_txs if exon_count(t, children) > 1]
    if multi_exon_txs and all_frameshifted:
        return 'pseudogene'

    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Flag pseudogenes in a gene GFF3.'
    )
    ap.add_argument('--gff3',                required=True)
    ap.add_argument('--repeats',             required=True, help='Repeat GFF3')
    ap.add_argument('--out',                 required=True)
    ap.add_argument('--max-repeat-coverage', type=float, default=0.80)
    args = ap.parse_args()

    headers, genes, transcripts, children, gene_order, tx_order = parse_gff3(args.gff3)
    repeat_map = parse_repeats(args.repeats)

    n_processed_pseudo = 0
    n_pseudo = 0

    for gid in gene_order:
        tx_ids = tx_order.get(gid, [])
        gene_rec = genes[gid]
        chrom = gene_rec.seqname
        repeat_ivs = repeat_map.get(chrom, [])

        biotype = classify_gene(
            gid, tx_ids, children, repeat_ivs, args.max_repeat_coverage
        )
        if biotype is None:
            continue

        # Update gene biotype attribute
        genes[gid].attrs = _set_attr(genes[gid].attrs, 'biotype', biotype)

        # Update transcript biotype attributes
        for txid in tx_ids:
            transcripts[txid].attrs = _set_attr(transcripts[txid].attrs, 'biotype', biotype)

        if biotype == 'processed_pseudogene':
            n_processed_pseudo += 1
        else:
            n_pseudo += 1

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
        f'detect_pseudogenes: {n_processed_pseudo} processed_pseudogenes, '
        f'{n_pseudo} pseudogenes flagged',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
