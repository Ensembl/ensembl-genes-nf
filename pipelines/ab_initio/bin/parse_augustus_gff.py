#!/usr/bin/env python3
"""
parse_augustus_gff.py — Convert Augustus GFF3 output to Ensembl-style GFF3.

Augustus outputs genes in a non-standard GFF3 format; this script normalises:
  - Extracts gene / mRNA / exon / CDS features
  - Assigns stable IDs (ab_initio_gene_NNNNNNNN, etc.)
  - Sets biotype=ab_initio
  - Filters by minimum gene length

Usage:
    parse_augustus_gff.py \\
        --input  augustus.gff \\
        --out    ab_initio.gff3 \\
        --min_gene_length 100
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class AugGene:
    seqname: str
    start:   int      # 1-based
    end:     int
    strand:  str
    gene_id: str
    transcripts: List['AugTranscript'] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class AugTranscript:
    seqname: str
    start:   int
    end:     int
    strand:  str
    tx_id:   str
    gene_id: str
    exons:   List[Tuple[int, int]] = field(default_factory=list)
    cdss:    List[Tuple[int, int]] = field(default_factory=list)


def _open(path: str):
    return gzip.open if path.endswith('.gz') else open


def parse_augustus_gff(path: str) -> List[AugGene]:
    """
    Parse Augustus GFF/GFF3 output.

    Augustus emits features in this order:
      gene   → ID=gX
      mRNA   → ID=gX.tY;Parent=gX   (or transcript in some modes)
      CDS    → Parent=gX.tY
      exon   → Parent=gX.tY

    Some Augustus versions emit # comment blocks before each gene.
    """
    genes:   Dict[str, AugGene]       = {}
    txs:     Dict[str, AugTranscript] = {}
    tx_exons: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    tx_cdss:  Dict[str, List[Tuple[int,int]]] = defaultdict(list)

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            seqname, _, feature, start, end, _, strand, _, attrs_raw = cols[:9]
            start, end = int(start), int(end)

            attrs = _parse_attrs(attrs_raw)

            if feature == 'gene':
                gid = attrs.get('ID') or attrs.get('gene_id', '')
                if gid:
                    genes[gid] = AugGene(
                        seqname=seqname, start=start, end=end,
                        strand=strand, gene_id=gid,
                    )

            elif feature in ('mRNA', 'transcript'):
                tid  = attrs.get('ID', '')
                gid  = attrs.get('Parent', attrs.get('gene_id', ''))
                if tid:
                    txs[tid] = AugTranscript(
                        seqname=seqname, start=start, end=end,
                        strand=strand, tx_id=tid, gene_id=gid,
                    )
                    if gid in genes:
                        genes[gid].transcripts.append(txs[tid])

            elif feature == 'exon':
                parent = attrs.get('Parent', '')
                if parent:
                    tx_exons[parent].append((start, end))

            elif feature == 'CDS':
                parent = attrs.get('Parent', '')
                if parent:
                    tx_cdss[parent].append((start, end))

    # Attach exons/CDS to transcripts
    for tid, tx in txs.items():
        tx.exons = sorted(tx_exons.get(tid, []))
        tx.cdss  = sorted(tx_cdss.get(tid, []))
        # Augustus sometimes lacks explicit exon records — derive from CDS
        if not tx.exons and tx.cdss:
            tx.exons = tx.cdss[:]

    return list(genes.values())


def _parse_attrs(attr_str: str) -> dict:
    """Parse GFF9 attribute field (key=value; or key "value";)."""
    attrs = {}
    # Try GFF3 style first (key=value)
    for m in re.finditer(r'([^=;\s]+)=([^;]*)', attr_str):
        attrs[m.group(1)] = m.group(2).strip()
    if not attrs:
        # GTF style (key "value")
        for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
            attrs[m.group(1)] = m.group(2)
    return attrs


def filter_genes(genes: List[AugGene], min_gene_length: int) -> List[AugGene]:
    return [g for g in genes if g.length >= min_gene_length]


def write_gff3(genes: List[AugGene], out_path: str) -> Tuple[int, int]:
    """Write Ensembl-style GFF3.  Returns (gene_count, transcript_count)."""
    gene_count = tx_count = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, gene in enumerate(genes, 1):
            gene_ensid = f'ab_initio_gene_{i:08d}'
            fh.write(
                f'{gene.seqname}\taugustus\tgene'
                f'\t{gene.start}\t{gene.end}'
                f'\t.\t{gene.strand}\t.'
                f'\tID={gene_ensid};Name={gene.gene_id};biotype=ab_initio\n'
            )
            gene_count += 1

            for j, tx in enumerate(gene.transcripts or [], 1):
                tx_ensid = f'{gene_ensid}_tx_{j:04d}'
                fh.write(
                    f'{tx.seqname}\taugustus\ttranscript'
                    f'\t{tx.start}\t{tx.end}'
                    f'\t.\t{tx.strand}\t.'
                    f'\tID={tx_ensid};Parent={gene_ensid};Name={tx.tx_id};'
                    f'biotype=ab_initio\n'
                )
                tx_count += 1

                for k, (es, ee) in enumerate(tx.exons, 1):
                    fh.write(
                        f'{tx.seqname}\taugustus\texon'
                        f'\t{es}\t{ee}'
                        f'\t.\t{tx.strand}\t.'
                        f'\tID={tx_ensid}_exon_{k};Parent={tx_ensid}\n'
                    )

                for k, (cs, ce) in enumerate(tx.cdss, 1):
                    fh.write(
                        f'{tx.seqname}\taugustus\tCDS'
                        f'\t{cs}\t{ce}'
                        f'\t.\t{tx.strand}\t0'
                        f'\tID={tx_ensid}_cds_{k};Parent={tx_ensid}\n'
                    )

    return gene_count, tx_count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',           required=True)
    ap.add_argument('--out',             required=True)
    ap.add_argument('--min_gene_length', type=int, default=100)
    args = ap.parse_args()

    genes    = parse_augustus_gff(args.input)
    filtered = filter_genes(genes, args.min_gene_length)
    g_count, tx_count = write_gff3(filtered, args.out)

    print(
        f'parse_augustus_gff: {len(genes)} raw genes; '
        f'{g_count} retained (min_len={args.min_gene_length}); '
        f'{tx_count} transcripts',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
