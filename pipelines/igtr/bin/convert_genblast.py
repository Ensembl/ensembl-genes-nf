#!/usr/bin/env python3
"""
convert_genblast.py — Convert GenBlast GFF output to Ensembl GFF3.

Reads GenBlast's genBlastG output (which has 'transcript' and 'coding_exon'
features), applies PID and coverage filters, assigns ig_gene/tr_gene biotypes
from the IMGT protein FASTA headers, and writes a proper GFF3 with
gene/transcript/exon hierarchy.

GenBlast GFF input format:
  seqname  genBlastG  transcript  start  end  score  strand  .
    ID=IGHV1-2*02-R1-1-A1;Name=IGHV1-2*02;PID=85.00;Coverage=98.36
  seqname  genBlastG  coding_exon  start  end  .  strand  0
    ID=IGHV1-2*02-R1-1-A1-E1;Parent=IGHV1-2*02-R1-1-A1

Output GFF3:
  seqname  genBlastG  gene        start  end  score  strand  .
    ID=igtr_gene_00000001;Name=IGHV1-2*02;biotype=ig_gene
  seqname  genBlastG  transcript  start  end  score  strand  .
    ID=igtr_transcript_00000001;Parent=igtr_gene_00000001;
    Name=IGHV1-2*02;PID=85.00;Coverage=98.36
  seqname  genBlastG  exon        start  end  .  strand  .
    ID=igtr_exon_00000001;Parent=igtr_transcript_00000001

Usage:
    convert_genblast.py --gff genblast.gff --proteins batch.fa
                        --out filtered.gff3
                        [--min-pid 70] [--min-coverage 80]
                        [--sample-id SAMPLE]
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

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
        """Parse a GFF/GFF3 line. Returns None for comments/blank lines."""
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            return None
        cols = line.split('\t')
        if len(cols) < 9:
            return None
        attrs = parse_attrs(cols[8])
        return cls(
            seqname=cols[0],
            source=cols[1],
            feature=cols[2],
            start=int(cols[3]),   # GFF is 1-based
            end=int(cols[4]),
            score=cols[5],
            strand=cols[6],
            frame=cols[7],
            attrs=attrs,
        )

    def to_gff3_line(self, attrs_str: str) -> str:
        return '\t'.join([
            self.seqname, self.source, self.feature,
            str(self.start), str(self.end),
            self.score, self.strand, self.frame,
            attrs_str,
        ])


def parse_attrs(attr_str: str) -> dict:
    """Parse GFF3 key=value; attribute string into dict."""
    attrs = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if '=' in part:
            k, _, v = part.partition('=')
            attrs[k.strip()] = v.strip()
    return attrs


def attrs_to_str(attrs: dict) -> str:
    return ';'.join(f'{k}={v}' for k, v in attrs.items())


# ---------------------------------------------------------------------------
# Biotype assignment from IMGT FASTA headers
# ---------------------------------------------------------------------------

def load_protein_biotypes(fasta_path: str) -> dict[str, str]:
    """
    Parse IMGT FASTA headers (format: >biotype|name|source) and return
    a mapping from protein name → ig_gene | tr_gene.

    Falls back to name-prefix heuristic if header doesn't match IMGT format.
    """
    biotype_map: dict[str, str] = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].rstrip()
            parts = header.split('|')
            if len(parts) >= 2:
                raw_biotype = parts[0].strip()
                name = parts[1].strip()
            else:
                # No pipe — use whole header as name, infer biotype from prefix
                name = parts[0].split()[0]
                raw_biotype = name
            biotype_map[name] = _infer_biotype(raw_biotype, name)
    return biotype_map


def _infer_biotype(raw: str, name: str) -> str:
    """Map IMGT biotype/name prefix to ig_gene or tr_gene."""
    s = raw.upper()
    if s.startswith('IG') or name.upper().startswith('IG'):
        return 'ig_gene'
    if s.startswith('TR') or name.upper().startswith('TR'):
        return 'tr_gene'
    return 'ig_gene'   # default to ig_gene for unknown


# ---------------------------------------------------------------------------
# GenBlast GFF parsing
# ---------------------------------------------------------------------------

def parse_genblast_gff(gff_path: str) -> tuple[dict, dict]:
    """
    Parse GenBlast GFF output.

    Returns:
        transcripts: { transcript_id: GffRecord }
        exons:       { transcript_id: [GffRecord, ...] }
    """
    transcripts: dict[str, GffRecord] = {}
    exons: dict[str, list] = {}

    with open(gff_path) as fh:
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
# Filtering
# ---------------------------------------------------------------------------

def passes_filter(rec: GffRecord, min_pid: float, min_cov: float) -> bool:
    """Return True if transcript passes PID and coverage thresholds."""
    try:
        pid = float(rec.attrs.get('PID', 0))
        cov = float(rec.attrs.get('Coverage', 0))
    except ValueError:
        return False
    return pid >= min_pid and cov >= min_cov


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_gff3(
    transcripts: dict,
    exons: dict,
    biotype_map: dict,
    out_path: str,
    sample_id: str,
    min_pid: float,
    min_cov: float,
) -> int:
    """
    Filter transcripts, assign biotypes, write GFF3 with gene/transcript/exon.
    Returns the number of transcript records written.
    """
    gene_counter = 0
    tx_counter = 0
    exon_counter = 0
    count = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')

        for tid, tx_rec in transcripts.items():
            if not passes_filter(tx_rec, min_pid, min_cov):
                continue

            tx_exons = sorted(exons.get(tid, []), key=lambda r: r.start)
            if not tx_exons:
                continue

            # Get protein name and biotype
            prot_name = tx_rec.attrs.get('Name', tid)
            biotype = biotype_map.get(prot_name, _infer_biotype('', prot_name))
            pid = tx_rec.attrs.get('PID', '0')
            cov = tx_rec.attrs.get('Coverage', '0')

            gene_counter += 1
            tx_counter += 1
            gene_id = f'{sample_id}_igtr_gene_{gene_counter:08d}'
            tx_id   = f'{sample_id}_igtr_transcript_{tx_counter:08d}'

            # Gene span covers all exons
            gene_start = tx_exons[0].start
            gene_end   = tx_exons[-1].end

            # Gene record
            gene_attrs = attrs_to_str({
                'ID':      gene_id,
                'Name':    prot_name,
                'biotype': biotype,
            })
            fh.write('\t'.join([
                tx_rec.seqname, 'genBlastG', 'gene',
                str(gene_start), str(gene_end),
                tx_rec.score, tx_rec.strand, '.',
                gene_attrs,
            ]) + '\n')

            # Transcript record
            tx_attrs = attrs_to_str({
                'ID':       tx_id,
                'Parent':   gene_id,
                'Name':     prot_name,
                'PID':      pid,
                'Coverage': cov,
                'biotype':  biotype,
            })
            fh.write('\t'.join([
                tx_rec.seqname, 'genBlastG', 'transcript',
                str(tx_rec.start), str(tx_rec.end),
                tx_rec.score, tx_rec.strand, '.',
                tx_attrs,
            ]) + '\n')

            # Exon records
            for exon_rec in tx_exons:
                exon_counter += 1
                exon_id = f'{sample_id}_igtr_exon_{exon_counter:08d}'
                exon_attrs = attrs_to_str({
                    'ID':     exon_id,
                    'Parent': tx_id,
                })
                fh.write('\t'.join([
                    exon_rec.seqname, 'genBlastG', 'exon',
                    str(exon_rec.start), str(exon_rec.end),
                    '.', exon_rec.strand, '.',
                    exon_attrs,
                ]) + '\n')

            count += 1

    return count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff',          required=True, help='GenBlast GFF output')
    parser.add_argument('--proteins',     required=True, help='Protein FASTA (IMGT format headers)')
    parser.add_argument('--out',          required=True, help='Output GFF3 path')
    parser.add_argument('--min-pid',      type=float, default=70.0, help='Min percent identity (default: 70)')
    parser.add_argument('--min-coverage', type=float, default=80.0, help='Min query coverage (default: 80)')
    parser.add_argument('--sample-id',    default='sample', help='Sample ID prefix for feature IDs')
    args = parser.parse_args()

    print(f'[convert_genblast] Loading biotypes from {args.proteins}', flush=True)
    biotype_map = load_protein_biotypes(args.proteins)
    print(f'[convert_genblast] Loaded {len(biotype_map)} protein biotypes', flush=True)

    print(f'[convert_genblast] Parsing GenBlast GFF: {args.gff}', flush=True)
    transcripts, exons = parse_genblast_gff(args.gff)
    print(f'[convert_genblast] Found {len(transcripts)} transcripts', flush=True)

    n = write_gff3(
        transcripts, exons, biotype_map,
        args.out, args.sample_id,
        args.min_pid, args.min_coverage,
    )
    print(f'[convert_genblast] Wrote {n} models to {args.out} '
          f'(PID≥{args.min_pid}%, Coverage≥{args.min_coverage}%)', flush=True)


if __name__ == '__main__':
    main()
