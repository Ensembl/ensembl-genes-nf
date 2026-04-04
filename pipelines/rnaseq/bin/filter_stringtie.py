#!/usr/bin/env python3
"""
filter_stringtie.py — Filter StringTie2 GTF output and convert to GFF3.

StringTie2 GTF columns:
  chr  source  feature  start  end  score  strand  frame  attributes

Features used: transcript (for model) + exon (child features)
Attributes: transcript_id, gene_id, cov (coverage), FPKM, TPM

Filtering thresholds:
  --min_coverage   : minimum transcript read coverage (cov attribute)
  --min_length     : minimum transcript length in bp
  --min_exons      : minimum number of exons

Output: Ensembl-style GFF3 with biotype=rnaseq_tissue or rnaseq_merged.

Usage:
    filter_stringtie.py \\
        --gtf      stringtie.gtf \\
        --out      rnaseq.gff3 \\
        --sample   SAMPLE001 \\
        --biotype  rnaseq_tissue \\
        --min_coverage 2.0 \\
        --min_length   200 \\
        --min_exons    1
"""

import argparse
import re
import sys
from dataclasses import dataclass, field
from typing import List, Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Transcript:
    seqname:     str
    start:       int
    end:         int
    strand:      str
    transcript_id: str
    gene_id:     str
    coverage:    float
    fpkm:        float
    tpm:         float
    exons:       List = field(default_factory=list)  # list of (start, end)

    @property
    def length(self) -> int:
        return sum(e - s + 1 for s, e in self.exons) if self.exons else self.end - self.start + 1

    @property
    def n_exons(self) -> int:
        return len(self.exons)


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_gtf_attrs(attr_str: str) -> dict:
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr_str)}


def _float_attr(attrs: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(attrs.get(key, default))
    except (ValueError, TypeError):
        return default


def _open(path: str):
    import gzip
    if path.endswith('.gz'):
        return gzip.open
    return open


def parse_stringtie_gtf(path: str) -> List[Transcript]:
    """Parse StringTie2 GTF and return list of Transcript objects."""
    txs = {}
    exons = {}

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ('transcript', 'exon'):
                continue

            attrs = _parse_gtf_attrs(cols[8])
            tx_id = attrs.get('transcript_id', '')
            if not tx_id:
                continue

            s, e = int(cols[3]), int(cols[4])

            if feature == 'transcript':
                txs[tx_id] = Transcript(
                    seqname      = cols[0],
                    start        = s,
                    end          = e,
                    strand       = cols[6],
                    transcript_id = tx_id,
                    gene_id      = attrs.get('gene_id', tx_id),
                    coverage     = _float_attr(attrs, 'cov'),
                    fpkm         = _float_attr(attrs, 'FPKM'),
                    tpm          = _float_attr(attrs, 'TPM'),
                )
                exons[tx_id] = []
            elif feature == 'exon' and tx_id in txs:
                exons[tx_id].append((s, e))

    # Attach exons
    for tx_id, tx in txs.items():
        tx.exons = sorted(exons.get(tx_id, []))

    return list(txs.values())


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_transcripts(
    txs:           List[Transcript],
    min_coverage:  float,
    min_length:    int,
    min_exons:     int,
) -> List[Transcript]:
    return [
        tx for tx in txs
        if tx.coverage  >= min_coverage
        and tx.length   >= min_length
        and tx.n_exons  >= min_exons
    ]


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_gff3(
    txs:      List[Transcript],
    out_path: str,
    sample:   str,
    biotype:  str,
) -> int:
    """Write gene/transcript/exon GFF3.  Returns number of transcripts written."""
    # Group by gene_id to write gene features
    from collections import defaultdict
    gene_txs = defaultdict(list)
    for tx in txs:
        gene_txs[tx.gene_id].append(tx)

    written = 0
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for gene_id, gene_transcripts in gene_txs.items():
            seqname = gene_transcripts[0].seqname
            strand  = gene_transcripts[0].strand
            g_start = min(tx.start for tx in gene_transcripts)
            g_end   = max(tx.end   for tx in gene_transcripts)
            gene_ensid = f'{sample}_rna_gene_{gene_id}'
            fh.write(
                f'{seqname}\tStringTie2\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.'
                f'\tID={gene_ensid};Name={gene_id};biotype={biotype}\n'
            )
            for tx in sorted(gene_transcripts, key=lambda t: t.start):
                tx_ensid = f'{sample}_rna_tx_{tx.transcript_id}'
                fh.write(
                    f'{tx.seqname}\tStringTie2\ttranscript\t{tx.start}\t{tx.end}'
                    f'\t{tx.coverage:.2f}\t{tx.strand}\t.'
                    f'\tID={tx_ensid};Parent={gene_ensid};Name={tx.transcript_id};'
                    f'biotype={biotype};cov={tx.coverage:.2f};FPKM={tx.fpkm:.4f};TPM={tx.tpm:.4f}\n'
                )
                for i, (es, ee) in enumerate(tx.exons, 1):
                    fh.write(
                        f'{tx.seqname}\tStringTie2\texon\t{es}\t{ee}'
                        f'\t.\t{tx.strand}\t.'
                        f'\tID={tx_ensid}_exon_{i};Parent={tx_ensid}\n'
                    )
                written += 1
    return written


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf',          required=True)
    ap.add_argument('--out',          required=True)
    ap.add_argument('--sample',       required=True)
    ap.add_argument('--biotype',      default='rnaseq_tissue')
    ap.add_argument('--min_coverage', type=float, default=2.0)
    ap.add_argument('--min_length',   type=int,   default=200)
    ap.add_argument('--min_exons',    type=int,   default=1)
    args = ap.parse_args()

    txs      = parse_stringtie_gtf(args.gtf)
    passing  = filter_transcripts(txs, args.min_coverage, args.min_length, args.min_exons)
    n        = write_gff3(passing, args.out, args.sample, args.biotype)

    print(f'filter_stringtie: {len(txs)} transcripts parsed, {n} passed filters',
          file=sys.stderr)


if __name__ == '__main__':
    main()
