#!/usr/bin/env python3
"""
classify_genblast.py — Convert GenBlast GFF to tiered Ensembl GFF3.

Reads GenBlast genBlastG output, applies PID/coverage filters, assigns
biotypes genblast_1..7 using the standard Ensembl classification tiers,
and writes GFF3 with gene/transcript/exon hierarchy.

Classification tiers (matching HiveClassifyTranscriptSupport 'standard'):
  1 — coverage>=95  AND identity>=90
  2 — coverage>=90  AND identity>=80
  3 — coverage>=90  AND identity>=60
  4 — coverage>=90  AND identity>=40
  5 — coverage>=80  AND identity>=20
  6 — coverage>=60  AND identity>=20
  7 — everything else (passes min thresholds but no tier match)

Usage:
    classify_genblast.py --gff genblast.gff --out classified.gff3
                         [--min-pid 30] [--min-coverage 50]
                         [--sample-id SAMPLE]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Classification tiers — (min_coverage, min_identity) → biotype
# ---------------------------------------------------------------------------

TIERS = [
    (95, 90, 'genblast_1'),
    (90, 80, 'genblast_2'),
    (90, 60, 'genblast_3'),
    (90, 40, 'genblast_4'),
    (80, 20, 'genblast_5'),
    (60, 20, 'genblast_6'),
]


def classify(pid: float, coverage: float) -> str:
    """Return genblast_N biotype for given PID and coverage values."""
    for min_cov, min_pid, biotype in TIERS:
        if coverage >= min_cov and pid >= min_pid:
            return biotype
    return 'genblast_7'


# ---------------------------------------------------------------------------
# GFF parsing (same pattern as convert_genblast.py)
# ---------------------------------------------------------------------------

def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if '=' in part:
            k, _, v = part.partition('=')
            attrs[k.strip()] = v.strip()
    return attrs


def attrs_to_str(attrs: dict) -> str:
    return ';'.join(f'{k}={v}' for k, v in attrs.items())


@dataclass
class GffRecord:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    attrs: dict = field(default_factory=dict)

    @classmethod
    def from_line(cls, line: str) -> Optional['GffRecord']:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            return None
        cols = line.split('\t')
        if len(cols) < 9:
            return None
        return cls(
            seqname=cols[0], source=cols[1], feature=cols[2],
            start=int(cols[3]), end=int(cols[4]),
            score=cols[5], strand=cols[6], frame=cols[7],
            attrs=parse_attrs(cols[8]),
        )


def parse_genblast_gff(path: str) -> tuple[dict, dict]:
    transcripts: dict = {}
    exons: dict = {}
    with open(path) as fh:
        for line in fh:
            rec = GffRecord.from_line(line)
            if rec is None:
                continue
            if rec.feature == 'transcript':
                tid = rec.attrs.get('ID', '')
                if tid:
                    transcripts[tid] = rec
            elif rec.feature == 'coding_exon':
                parent = rec.attrs.get('Parent', '')
                if parent:
                    exons.setdefault(parent, []).append(rec)
    return transcripts, exons


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_classified_gff3(
    transcripts: dict,
    exons: dict,
    out_path: str,
    sample_id: str,
    min_pid: float,
    min_cov: float,
) -> dict:
    """Write GFF3 with genblast_N biotypes. Returns count per tier."""
    counts: dict = {}
    gene_n = tx_n = exon_n = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for tid, tx in transcripts.items():
            try:
                pid = float(tx.attrs.get('PID', 0))
                cov = float(tx.attrs.get('Coverage', 0))
            except ValueError:
                continue
            if pid < min_pid or cov < min_cov:
                continue

            tx_exons = sorted(exons.get(tid, []), key=lambda r: r.start)
            if not tx_exons:
                continue

            biotype = classify(pid, cov)
            counts[biotype] = counts.get(biotype, 0) + 1

            prot_name = tx.attrs.get('Name', tid)
            gene_n += 1; tx_n += 1
            gene_id = f'{sample_id}_gbh_gene_{gene_n:08d}'
            tx_id   = f'{sample_id}_gbh_transcript_{tx_n:08d}'

            gene_start = tx_exons[0].start
            gene_end   = tx_exons[-1].end

            fh.write('\t'.join([
                tx.seqname, 'genBlastG', 'gene',
                str(gene_start), str(gene_end),
                tx.score, tx.strand, '.',
                attrs_to_str({'ID': gene_id, 'Name': prot_name, 'biotype': biotype}),
            ]) + '\n')

            fh.write('\t'.join([
                tx.seqname, 'genBlastG', 'transcript',
                str(tx.start), str(tx.end),
                tx.score, tx.strand, '.',
                attrs_to_str({
                    'ID': tx_id, 'Parent': gene_id, 'Name': prot_name,
                    'PID': str(pid), 'Coverage': str(cov), 'biotype': biotype,
                }),
            ]) + '\n')

            for exon in tx_exons:
                exon_n += 1
                fh.write('\t'.join([
                    exon.seqname, 'genBlastG', 'exon',
                    str(exon.start), str(exon.end),
                    '.', exon.strand, '.',
                    attrs_to_str({'ID': f'{sample_id}_gbh_exon_{exon_n:08d}', 'Parent': tx_id}),
                ]) + '\n')

    return counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff',          required=True)
    parser.add_argument('--out',          required=True)
    parser.add_argument('--min-pid',      type=float, default=30.0)
    parser.add_argument('--min-coverage', type=float, default=50.0)
    parser.add_argument('--sample-id',    default='sample')
    args = parser.parse_args()

    txs, exons = parse_genblast_gff(args.gff)
    print(f'[classify_genblast] {len(txs)} transcripts parsed from {args.gff}', flush=True)

    counts = write_classified_gff3(
        txs, exons, args.out, args.sample_id,
        args.min_pid, args.min_coverage,
    )
    total = sum(counts.values())
    print(f'[classify_genblast] Wrote {total} models to {args.out}', flush=True)
    for bt in sorted(counts):
        print(f'  {bt}: {counts[bt]}', flush=True)


if __name__ == '__main__':
    main()
