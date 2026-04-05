#!/usr/bin/env python3
"""
flag_selenoproteins.py — Flag genes matching a selenoprotein database.

When the selenoprotein FASTA is the NO_FILE sentinel (name == 'NO_FILE'),
the input GFF3 is copied to the output unchanged.

Otherwise, gene Name attributes are matched against sequence IDs in the FASTA.
Matching genes have their biotype set to 'selenoprotein'.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set


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


def parse_selenoprotein_names(fasta_path: str) -> Set[str]:
    """
    Return the set of sequence IDs from a FASTA file.
    The ID is the first whitespace-delimited token of the defline (without '>').
    """
    names: Set[str] = set()
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                names.add(seq_id)
    return names


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Flag selenoprotein genes in a GFF3.'
    )
    ap.add_argument('--gff3',     required=True, help='Input GFF3')
    ap.add_argument('--proteins', required=True, help='Selenoprotein FASTA or NO_FILE sentinel')
    ap.add_argument('--out',      required=True, help='Output GFF3')
    args = ap.parse_args()

    # Check for NO_FILE sentinel
    if os.path.basename(args.proteins) == 'NO_FILE':
        shutil.copy(args.gff3, args.out)
        print('flag_selenoproteins: no selenoprotein FASTA provided — input copied unchanged',
              file=sys.stderr)
        return

    seleno_names = parse_selenoprotein_names(args.proteins)
    headers, genes, transcripts, children, gene_order, tx_order = parse_gff3(args.gff3)

    n_flagged = 0

    for gid in gene_order:
        gene_rec = genes[gid]
        attrs = _parse_attrs(gene_rec.attrs)
        gene_name = attrs.get('Name', attrs.get('ID', ''))

        if gene_name in seleno_names:
            genes[gid].attrs = _set_attr(genes[gid].attrs, 'biotype', 'selenoprotein')
            for txid in tx_order.get(gid, []):
                transcripts[txid].attrs = _set_attr(
                    transcripts[txid].attrs, 'biotype', 'selenoprotein'
                )
            n_flagged += 1

    with open(args.out, 'w') as fh:
        for h in headers:
            fh.write(h)
        for gid in gene_order:
            fh.write(_rebuild_raw(genes[gid]) + '\n')
            for txid in tx_order.get(gid, []):
                fh.write(_rebuild_raw(transcripts[txid]) + '\n')
                for child in children.get(txid, []):
                    fh.write(_rebuild_raw(child) + '\n')

    print(f'flag_selenoproteins: {n_flagged} selenoprotein genes flagged', file=sys.stderr)


if __name__ == '__main__':
    main()
