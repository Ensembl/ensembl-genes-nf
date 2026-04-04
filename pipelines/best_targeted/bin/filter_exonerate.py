#!/usr/bin/env python3
"""
filter_exonerate.py — Parse exonerate GFF2 output and filter by coverage/percent-id.

Exonerate (--showtargetgff) produces GFF2-style output.  We extract gene-level
records and the supporting 'exon' features, compute alignment coverage and pid
from the exonerate-specific tags in the score column, then write GFF3.

Usage:
    filter_exonerate.py \\
        --exonerate_gff <file.gff> \\
        --out            <file.gff3> \\
        --query_type     cdna|protein \\
        --min_coverage   50 \\
        --min_pid        50 \\
        --sample_id      CHUNK01
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
class ExonerateHit:
    seqname:   str
    start:     int        # 1-based, always start < end
    end:       int
    strand:    str
    score:     float
    coverage:  float
    pid:       float
    query_id:  str
    biotype:   str
    exons:     List       = field(default_factory=list)  # list of (start, end, strand)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_GENE_RE = re.compile(
    r'^(\S+)\s+exonerate:(\S+)\s+gene\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+\.\s+(.+)$'
)
_EXON_RE = re.compile(
    r'^(\S+)\s+exonerate:\S+\s+exon\s+(\d+)\s+(\d+)\s+\S+\s+([+-])\s+\.\s+(.+)$'
)
_ATTR_RE = re.compile(r'(\w+)\s+([^;]+)')


def _parse_attrs(attr_str: str) -> dict:
    return {m.group(1): m.group(2).strip() for m in _ATTR_RE.finditer(attr_str)}


def _gff2_coords(start: str, end: str):
    """Return 1-based (start, end) with start <= end."""
    s, e = int(start), int(end)
    if s > e:
        s, e = e, s
    return s, e


def parse_exonerate_gff(path: str, query_type: str) -> List[ExonerateHit]:
    """
    Parse exonerate --showtargetgff output.

    Exonerate GFF2 uses:
      score column  → alignment score (integer)
      attribute 'Query' → query sequence name
      'percentidentity', 'queryCoverage' may appear in vulgar/cigar lines but
      NOT in the --showtargetgff output.  We therefore compute coverage from
      the sum of aligned exon lengths vs the query length encoded in the
      sequence header (if available), and pid from the score heuristic.

    For cdna2genome:  exonerate emits  'gene', 'cds', 'exon', 'intron', 'utr' etc.
    For protein2genome: 'gene', 'similarity', 'exon'.

    We only need 'gene' + 'exon' features here.
    """
    hits: List[ExonerateHit] = []
    current: Optional[ExonerateHit] = None

    biotype = 'cdna_alignment' if query_type == 'cdna' else 'protein_alignment'

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) < 9:
                continue

            feature = cols[2]

            if feature == 'gene':
                m = _GENE_RE.match(line)
                if not m:
                    continue
                seqname = m.group(1)
                model   = m.group(2)    # e.g. cdna2genome
                s, e    = _gff2_coords(m.group(3), m.group(4))
                score   = float(m.group(5)) if m.group(5) != '.' else 0.0
                strand  = m.group(6)
                attrs   = _parse_attrs(m.group(7))
                query_id = attrs.get('Query', attrs.get('sequence', 'unknown'))
                # strip any trailing coordinate suffix exonerate sometimes adds
                query_id = query_id.split()[0]

                current = ExonerateHit(
                    seqname=seqname, start=s, end=e, strand=strand,
                    score=score, coverage=0.0, pid=0.0,
                    query_id=query_id, biotype=biotype,
                )
                hits.append(current)

            elif feature == 'exon' and current is not None:
                m = _EXON_RE.match(line)
                if m:
                    s, e = _gff2_coords(m.group(2), m.group(3))
                    strand = m.group(4)
                    current.exons.append((s, e, strand))

    return hits


def _open(path: str):
    import gzip
    if path.endswith('.gz'):
        return gzip.open
    return open


# ---------------------------------------------------------------------------
# Coverage / PID estimation
# ---------------------------------------------------------------------------

def _compute_stats(hit: ExonerateHit, query_lengths: dict) -> None:
    """
    Fill hit.coverage and hit.pid in-place.

    Exonerate --showtargetgff doesn't embed pid/coverage in GFF attributes,
    so we derive approximate values:
      coverage = sum(exon_lengths) / query_length  (if query length known)
      pid      = approximated from alignment score per aligned base
                 (score / aligned_bases * some_factor)

    When query lengths aren't known we fall back to genomic span coverage.
    This is intentionally approximate — downstream best-pick selection uses
    these values for ranking; the absolute numbers are less critical.
    """
    aligned_bases = sum(e - s + 1 for s, e, _ in hit.exons) if hit.exons else (hit.end - hit.start + 1)

    qlen = query_lengths.get(hit.query_id, 0)
    if qlen > 0:
        hit.coverage = min(100.0, aligned_bases / qlen * 100.0)
    else:
        # Fall back: use genomic span (underestimate for intron-containing genes)
        span = hit.end - hit.start + 1
        hit.coverage = min(100.0, aligned_bases / span * 100.0) if span > 0 else 0.0

    # Approximate pid from score: exonerate cdna2genome gives ~5 pts per match,
    # protein2genome ~10 pts per match.  We normalise to yield a rough pid.
    pts_per_match = 10.0 if 'protein' in hit.biotype else 5.0
    if aligned_bases > 0 and pts_per_match > 0:
        hit.pid = min(100.0, hit.score / (aligned_bases * pts_per_match) * 100.0)
    else:
        hit.pid = 0.0


def load_query_lengths(fasta_path: Optional[str]) -> dict:
    """Return {seq_id: length} from a FASTA file (supports .gz)."""
    if not fasta_path:
        return {}
    lengths = {}
    opener = _open(fasta_path)
    current_id = None
    current_len = 0
    with opener(fasta_path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            elif current_id:
                current_len += len(line)
    if current_id:
        lengths[current_id] = current_len
    return lengths


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_gff3(hits: List[ExonerateHit], out_path: str, sample_id: str) -> int:
    """Write GFF3 gene/transcript/exon hierarchy.  Returns number of models written."""
    written = 0
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, hit in enumerate(hits, 1):
            gene_id = f'{sample_id}_bt_gene_{i:08d}'
            tx_id   = f'{sample_id}_bt_tx_{i:08d}'
            attrs_gene = (f'ID={gene_id};Name={hit.query_id};'
                          f'biotype={hit.biotype};'
                          f'score={hit.score:.1f};'
                          f'coverage={hit.coverage:.1f};'
                          f'pid={hit.pid:.1f}')
            fh.write(f'{hit.seqname}\texonerate\tgene\t{hit.start}\t{hit.end}\t'
                     f'{hit.score:.1f}\t{hit.strand}\t.\t{attrs_gene}\n')
            attrs_tx = (f'ID={tx_id};Parent={gene_id};Name={hit.query_id};'
                        f'biotype={hit.biotype}')
            fh.write(f'{hit.seqname}\texonerate\ttranscript\t{hit.start}\t{hit.end}\t'
                     f'{hit.score:.1f}\t{hit.strand}\t.\t{attrs_tx}\n')
            for j, (es, ee, estrand) in enumerate(hit.exons, 1):
                exon_id = f'{tx_id}_exon_{j}'
                fh.write(f'{hit.seqname}\texonerate\texon\t{es}\t{ee}\t.\t'
                         f'{estrand}\t.\tID={exon_id};Parent={tx_id}\n')
            written += 1
    return written


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exonerate_gff', required=True)
    ap.add_argument('--out',           required=True)
    ap.add_argument('--query_type',    required=True, choices=['cdna', 'protein'])
    ap.add_argument('--min_coverage',  type=float, default=50.0)
    ap.add_argument('--min_pid',       type=float, default=50.0)
    ap.add_argument('--query_fasta',   default=None,
                    help='FASTA of query sequences (for accurate coverage)')
    ap.add_argument('--sample_id',     default='chunk')
    args = ap.parse_args()

    query_lengths = load_query_lengths(args.query_fasta)
    hits = parse_exonerate_gff(args.exonerate_gff, args.query_type)

    # compute stats and filter
    for hit in hits:
        _compute_stats(hit, query_lengths)

    passing = [h for h in hits
               if h.coverage >= args.min_coverage and h.pid >= args.min_pid]

    n = write_gff3(passing, args.out, args.sample_id)
    print(f'filter_exonerate: {len(hits)} hits parsed, {n} passed filters '
          f'(cov>={args.min_coverage}, pid>={args.min_pid})', file=sys.stderr)


if __name__ == '__main__':
    main()
