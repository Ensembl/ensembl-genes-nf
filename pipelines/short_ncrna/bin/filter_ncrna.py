#!/usr/bin/env python3
"""
filter_ncrna.py — Filter cmsearch tblout and BLASTN tabular output for ncRNA annotation.

Reads cmsearch --tblout output (Rfam search) and optionally a BLASTN tabular
output (miRBase search), applies thresholds, and writes a filtered GFF3 with
proper ncRNA biotypes.

Biotype assignment from Rfam CM names (simplified mapping):
  miRNA       — RF00001..miRNAs
  snRNA       — snRNA; families
  snoRNA/scaRNA — snoRNA; families
  rRNA        — rRNA; families (LSU, SSU, 5S, 5.8S)
  tRNA        — tRNA families
  ribozyme    — ribozyme; families
  misc_RNA    — everything else

cmsearch tblout columns (space-delimited, comments start with #):
  target_name, target_accession, query_name, query_accession, mdl, mdl_from,
  mdl_to, seq_from, seq_to, strand, trunc, pass, gc, bias, score, e_value,
  inc, description...

BLASTN tabular (-outfmt 6):
  qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
  sstart, send, evalue, bitscore

Usage:
    filter_ncrna.py --tblout rfam.tblout [--blast mirna.blast]
                    --out ncrna.gff3
                    [--min-score 0] [--max-evalue 0.01]
                    [--blast-min-pid 80] [--blast-max-evalue 0.01]
                    [--sample-id SAMPLE]
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Rfam biotype lookup — CM name prefix → biotype
# ---------------------------------------------------------------------------

# Ordered by specificity (scaRNA before snoRNA before snRNA)
RFAM_BIOTYPE_PATTERNS = [
    (re.compile(r'scaRNA',    re.I), 'scaRNA'),
    (re.compile(r'snoRNA|SNORD|SNORA|sno_', re.I), 'snoRNA'),
    (re.compile(r'snRNA|U1$|U2$|U4$|U5$|U6', re.I), 'snRNA'),
    (re.compile(r'miRNA|mir[-_]', re.I), 'miRNA'),
    (re.compile(r'rRNA|LSU_|SSU_|5S_|5_8S', re.I), 'rRNA'),
    (re.compile(r'tRNA',      re.I), 'tRNA'),
    (re.compile(r'Vault',     re.I), 'vault_RNA'),
    (re.compile(r'Y_RNA',     re.I), 'Y_RNA'),
    (re.compile(r'RNase|RNaseP', re.I), 'RNase_MRP_RNA'),
    (re.compile(r'ribozyme|hammerhead|HDV|GIR', re.I), 'ribozyme'),
    (re.compile(r'SRP',       re.I), 'SRP_RNA'),
    (re.compile(r'telomerase|TERC', re.I), 'telomerase_RNA'),
]


def rfam_biotype(cm_name: str) -> str:
    """Map Rfam CM name to Ensembl biotype string."""
    for pattern, biotype in RFAM_BIOTYPE_PATTERNS:
        if pattern.search(cm_name):
            return biotype
    return 'misc_RNA'


# ---------------------------------------------------------------------------
# cmsearch tblout parsing
# ---------------------------------------------------------------------------

@dataclass
class CmsearchHit:
    seqname: str
    cm_name: str
    cm_acc: str
    seq_from: int
    seq_to: int
    strand: str
    score: float
    evalue: float
    biotype: str = ''

    @classmethod
    def from_line(cls, line: str) -> Optional['CmsearchHit']:
        """Parse one non-comment line from cmsearch --tblout output."""
        if line.startswith('#') or not line.strip():
            return None
        # tblout is whitespace-delimited; first 18 fields matter
        cols = line.split()
        if len(cols) < 17:
            return None
        try:
            # Fields: target(0) tacc(1) query(2) qacc(3) mdl(4)
            #         mdl_from(5) mdl_to(6) seq_from(7) seq_to(8)
            #         strand(9) trunc(10) pass(11) gc(12) bias(13)
            #         score(14) evalue(15) inc(16)
            seq_from = int(cols[7])
            seq_to   = int(cols[8])
            strand   = cols[9]
            score    = float(cols[14])
            evalue   = float(cols[15])
            cm_name  = cols[2]
            cm_acc   = cols[3]
            seqname  = cols[0]
        except (ValueError, IndexError):
            return None

        # Normalise strand sign and coordinates
        if strand == '-' and seq_from > seq_to:
            seq_from, seq_to = seq_to, seq_from

        return cls(
            seqname=seqname,
            cm_name=cm_name,
            cm_acc=cm_acc,
            seq_from=seq_from,
            seq_to=seq_to,
            strand=strand,
            score=score,
            evalue=evalue,
            biotype=rfam_biotype(cm_name),
        )


def parse_tblout(path: str, max_evalue: float, min_score: float) -> list[CmsearchHit]:
    hits: list[CmsearchHit] = []
    # Support gzipped tblout (cmsearch can output .tbl.gz)
    opener = _open(path)
    with opener(path) as fh:
        for line in fh:
            if isinstance(line, bytes):
                line = line.decode()
            hit = CmsearchHit.from_line(line)
            if hit and hit.evalue <= max_evalue and hit.score >= min_score:
                hits.append(hit)
    return hits


# ---------------------------------------------------------------------------
# BLASTN tabular parsing (miRBase miRNA search)
# ---------------------------------------------------------------------------

@dataclass
class BlastHit:
    qseqid: str    # miRNA name
    sseqid: str    # genome sequence
    pident: float
    sstart: int
    send: int
    evalue: float
    bitscore: float
    strand: str = '+'  # derived from sstart/send


def parse_blast_tabular(path: str, min_pident: float, max_evalue: float) -> list[BlastHit]:
    hits: list[BlastHit] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 12:
                continue
            try:
                pident   = float(cols[2])
                sstart   = int(cols[8])
                send     = int(cols[9])
                evalue   = float(cols[10])
                bitscore = float(cols[11])
            except ValueError:
                continue
            if pident < min_pident or evalue > max_evalue:
                continue
            strand = '+' if sstart <= send else '-'
            if strand == '-':
                sstart, send = send, sstart
            hits.append(BlastHit(
                qseqid=cols[0], sseqid=cols[1],
                pident=pident, sstart=sstart, send=send,
                evalue=evalue, bitscore=bitscore, strand=strand,
            ))
    return hits


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_gff3(
    rfam_hits: list[CmsearchHit],
    blast_hits: list[BlastHit],
    out_path: str,
    sample_id: str,
) -> dict:
    gene_n = tx_n = 0
    counts: dict = {}

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')

        for hit in rfam_hits:
            gene_n += 1; tx_n += 1
            gene_id = f'{sample_id}_ncrna_gene_{gene_n:08d}'
            tx_id   = f'{sample_id}_ncrna_transcript_{tx_n:08d}'
            bt = hit.biotype
            counts[bt] = counts.get(bt, 0) + 1

            fh.write('\t'.join([
                hit.seqname, 'cmsearch', 'gene',
                str(hit.seq_from), str(hit.seq_to),
                str(hit.score), hit.strand, '.',
                f'ID={gene_id};Name={hit.cm_name};biotype={bt};rfam_acc={hit.cm_acc}',
            ]) + '\n')
            fh.write('\t'.join([
                hit.seqname, 'cmsearch', 'transcript',
                str(hit.seq_from), str(hit.seq_to),
                str(hit.score), hit.strand, '.',
                f'ID={tx_id};Parent={gene_id};Name={hit.cm_name};biotype={bt};'
                f'evalue={hit.evalue};score={hit.score}',
            ]) + '\n')
            fh.write('\t'.join([
                hit.seqname, 'cmsearch', 'exon',
                str(hit.seq_from), str(hit.seq_to),
                '.', hit.strand, '.',
                f'ID={sample_id}_ncrna_exon_{gene_n:08d};Parent={tx_id}',
            ]) + '\n')

        for hit in blast_hits:
            gene_n += 1; tx_n += 1
            gene_id = f'{sample_id}_ncrna_gene_{gene_n:08d}'
            tx_id   = f'{sample_id}_ncrna_transcript_{tx_n:08d}'
            bt = 'miRNA'
            counts[bt] = counts.get(bt, 0) + 1

            fh.write('\t'.join([
                hit.sseqid, 'blastn', 'gene',
                str(hit.sstart), str(hit.send),
                str(hit.bitscore), hit.strand, '.',
                f'ID={gene_id};Name={hit.qseqid};biotype={bt}',
            ]) + '\n')
            fh.write('\t'.join([
                hit.sseqid, 'blastn', 'transcript',
                str(hit.sstart), str(hit.send),
                str(hit.bitscore), hit.strand, '.',
                f'ID={tx_id};Parent={gene_id};Name={hit.qseqid};biotype={bt};'
                f'pident={hit.pident};evalue={hit.evalue}',
            ]) + '\n')
            fh.write('\t'.join([
                hit.sseqid, 'blastn', 'exon',
                str(hit.sstart), str(hit.send),
                '.', hit.strand, '.',
                f'ID={sample_id}_ncrna_exon_{gene_n:08d};Parent={tx_id}',
            ]) + '\n')

    return counts


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open(path: str):
    import gzip
    if path.endswith('.gz'):
        return gzip.open
    return open


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tblout',          required=True, help='cmsearch --tblout file (plain or .gz)')
    parser.add_argument('--blast',           default=None,  help='BLASTN tabular (-outfmt 6) for miRBase')
    parser.add_argument('--out',             required=True)
    parser.add_argument('--min-score',       type=float, default=0.0,  help='Min cmsearch bit score')
    parser.add_argument('--max-evalue',      type=float, default=0.01, help='Max cmsearch e-value')
    parser.add_argument('--blast-min-pid',   type=float, default=80.0, help='Min BLASTN identity for miRNA')
    parser.add_argument('--blast-max-evalue',type=float, default=0.01, help='Max BLASTN e-value for miRNA')
    parser.add_argument('--sample-id',       default='sample')
    args = parser.parse_args()

    rfam_hits = parse_tblout(args.tblout, args.max_evalue, args.min_score)
    print(f'[filter_ncrna] {len(rfam_hits)} Rfam hits after filtering', flush=True)

    blast_hits: list[BlastHit] = []
    if args.blast:
        blast_hits = parse_blast_tabular(args.blast, args.blast_min_pid, args.blast_max_evalue)
        print(f'[filter_ncrna] {len(blast_hits)} miRNA BLAST hits after filtering', flush=True)

    counts = write_gff3(rfam_hits, blast_hits, args.out, args.sample_id)
    total = sum(counts.values())
    print(f'[filter_ncrna] Wrote {total} ncRNA models to {args.out}', flush=True)
    for bt in sorted(counts):
        print(f'  {bt}: {counts[bt]}', flush=True)


if __name__ == '__main__':
    main()
