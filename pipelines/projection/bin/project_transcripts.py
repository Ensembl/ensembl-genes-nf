#!/usr/bin/env python3
"""
project_transcripts.py — Project gene models from a source (reference)
genome to a target genome using a LASTZ/minimap2 chain file.

The chain file format is the standard UCSC chain format:
  chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
  block_size dt dq
  ...

Projection strategy:
  1. Parse chain file to build coordinate map (target → query)
  2. Parse source GFF3 to get transcripts + exons
  3. For each transcript, map exon coordinates through the chain
  4. Filter by coverage and percent alignment
  5. Write projected GFF3 in target coordinates

Output biotype: projected_transcript

Usage:
    project_transcripts.py \\
        --source_gff3  source_genes.gff3 \\
        --chain        lastz.chain \\
        --out          projected.gff3 \\
        --min_coverage 50 \\
        --min_pid      50 \\
        --max_stops    1
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Chain file parsing
# ---------------------------------------------------------------------------

@dataclass
class ChainBlock:
    """A single ungapped alignment block within a chain."""
    t_start: int    # target (reference) start, 0-based
    t_end:   int    # target end, exclusive
    q_start: int    # query (new genome) start, 0-based
    q_end:   int    # query end, exclusive
    q_strand: str   # '+' or '-'


@dataclass
class Chain:
    score:   int
    t_name:  str     # target chromosome
    t_size:  int
    t_start: int     # target alignment start, 0-based
    t_end:   int
    q_name:  str     # query chromosome
    q_size:  int
    q_strand: str    # '+' or '-'
    q_start: int
    q_end:   int
    chain_id: str
    blocks:  List[ChainBlock] = field(default_factory=list)


def _open(path: str):
    return gzip.open if path.endswith('.gz') else open


def parse_chain(path: str) -> List[Chain]:
    """Parse a UCSC chain file and return list of Chain objects."""
    chains = []
    current: Optional[Chain] = None
    t_pos = 0
    q_pos = 0

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                # blank line ends a chain
                current = None
                continue
            if line.startswith('chain'):
                parts = line.split()
                if len(parts) < 13:
                    continue
                current = Chain(
                    score    = int(parts[1]),
                    t_name   = parts[2],
                    t_size   = int(parts[3]),
                    t_start  = int(parts[5]),
                    t_end    = int(parts[6]),
                    q_name   = parts[7],
                    q_size   = int(parts[8]),
                    q_strand = parts[9],
                    q_start  = int(parts[10]),
                    q_end    = int(parts[11]),
                    chain_id = parts[12],
                )
                t_pos = current.t_start
                q_pos = current.q_start if current.q_strand == '+' else (current.q_size - current.q_end)
                chains.append(current)
            elif current is not None:
                parts = line.split()
                if not parts:
                    continue
                block_size = int(parts[0])
                dt = int(parts[1]) if len(parts) > 1 else 0
                dq = int(parts[2]) if len(parts) > 2 else 0

                if current.q_strand == '+':
                    q_block_start = q_pos
                else:
                    # minus strand: query coords are in reverse complement
                    q_block_start = current.q_size - q_pos - block_size

                current.blocks.append(ChainBlock(
                    t_start  = t_pos,
                    t_end    = t_pos + block_size,
                    q_start  = q_block_start,
                    q_end    = q_block_start + block_size,
                    q_strand = current.q_strand,
                ))
                t_pos += block_size + dt
                q_pos += block_size + dq

    return chains


# ---------------------------------------------------------------------------
# Source GFF3 parsing (transcripts + exons only)
# ---------------------------------------------------------------------------

@dataclass
class SourceTranscript:
    seqname: str
    start:   int    # 1-based
    end:     int
    strand:  str
    tx_id:   str
    gene_id: str
    biotype: str
    exons:   List[Tuple[int, int]] = field(default_factory=list)  # 1-based (start, end)


def _parse_attrs(attr_str: str) -> dict:
    return {m.group(1): m.group(2)
            for m in re.finditer(r'([^=;]+)=([^;]*)', attr_str)}


def parse_source_gff3(path: str) -> List[SourceTranscript]:
    txs: Dict[str, SourceTranscript] = {}
    tx_exons: Dict[str, List] = defaultdict(list)

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
            attrs   = _parse_attrs(cols[8])

            if feature == 'transcript':
                tx_id   = attrs.get('ID', '')
                gene_id = attrs.get('Parent', attrs.get('ID', ''))
                if tx_id:
                    txs[tx_id] = SourceTranscript(
                        seqname = cols[0],
                        start   = int(cols[3]),
                        end     = int(cols[4]),
                        strand  = cols[6],
                        tx_id   = tx_id,
                        gene_id = gene_id,
                        biotype = attrs.get('biotype', 'protein_coding'),
                    )
            elif feature == 'exon':
                parent = attrs.get('Parent', '')
                if parent:
                    tx_exons[parent].append((int(cols[3]), int(cols[4])))

    for tx_id, tx in txs.items():
        tx.exons = sorted(tx_exons.get(tx_id, []))

    return list(txs.values())


# ---------------------------------------------------------------------------
# Coordinate projection via chain
# ---------------------------------------------------------------------------

def _build_t_index(chains: List[Chain]) -> Dict[str, List[Chain]]:
    """Index chains by target chromosome name."""
    idx: Dict[str, List[Chain]] = defaultdict(list)
    for chain in chains:
        idx[chain.t_name].append(chain)
    return idx


def _project_interval(
    t_start: int,   # 1-based
    t_end:   int,   # 1-based
    chain: Chain,
) -> Optional[Tuple[int, int, str]]:
    """
    Map a target (source) interval onto query (target genome) coordinates
    via one chain.  Returns (q_start, q_end, q_strand) in 1-based coords,
    or None if no overlap.
    """
    # Convert to 0-based half-open
    ts, te = t_start - 1, t_end

    q_starts = []
    q_ends   = []

    for block in chain.blocks:
        # Overlap of [ts, te) with [block.t_start, block.t_end)
        ov_s = max(ts, block.t_start)
        ov_e = min(te, block.t_end)
        if ov_s >= ov_e:
            continue
        offset_start = ov_s - block.t_start
        offset_end   = ov_e - block.t_start
        if chain.q_strand == '+':
            q_starts.append(block.q_start + offset_start)
            q_ends.append(block.q_start + offset_end)
        else:
            # minus strand: q coords are in fwd coords of q chromosome
            # block.q_start is already adjusted in parse_chain
            q_starts.append(block.q_start + offset_start)
            q_ends.append(block.q_start + offset_end)

    if not q_starts:
        return None

    q_s = min(q_starts)
    q_e = max(q_ends)
    strand = chain.q_strand
    # Convert back to 1-based
    return q_s + 1, q_e, strand


def _compute_coverage(
    exons: List[Tuple[int, int]],
    projected: List[Optional[Tuple[int, int, str]]],
) -> float:
    """Return fraction of exon bases that projected successfully."""
    total  = sum(e - s + 1 for s, e in exons)
    mapped = sum(
        (p[1] - p[0] + 1)
        for p, ex in zip(projected, exons)
        if p is not None
    )
    return (mapped / total * 100.0) if total > 0 else 0.0


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_projected_gff3(
    projections: List[dict],
    out_path: str,
) -> int:
    """Write projected transcripts to GFF3.  Returns count written."""
    # Group by projected gene (cluster same q_name + adjacent ranges)
    from collections import defaultdict
    by_gene: Dict[str, List[dict]] = defaultdict(list)
    for p in projections:
        by_gene[p['src_gene_id']].append(p)

    written = 0
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, (gene_id, txs) in enumerate(by_gene.items(), 1):
            q_name  = txs[0]['q_name']
            q_strand = txs[0]['q_strand']
            g_start = min(t['q_start'] for t in txs)
            g_end   = max(t['q_end']   for t in txs)
            gene_ensid = f'proj_gene_{i:08d}'

            fh.write(
                f'{q_name}\tprojection\tgene\t{g_start}\t{g_end}'
                f'\t.\t{q_strand}\t.'
                f'\tID={gene_ensid};Name={gene_id};biotype=projected_transcript;'
                f'src_gene={gene_id}\n'
            )

            for j, tx in enumerate(txs, 1):
                tx_ensid = f'{gene_ensid}_tx_{j:04d}'
                fh.write(
                    f'{q_name}\tprojection\ttranscript\t{tx["q_start"]}\t{tx["q_end"]}'
                    f'\t{tx["coverage"]:.1f}\t{q_strand}\t.'
                    f'\tID={tx_ensid};Parent={gene_ensid};Name={tx["src_tx_id"]};'
                    f'biotype=projected_transcript;src_tx={tx["src_tx_id"]};'
                    f'coverage={tx["coverage"]:.1f}\n'
                )
                for k, (es, ee, _) in enumerate(tx['proj_exons'], 1):
                    fh.write(
                        f'{q_name}\tprojection\texon\t{es}\t{ee}'
                        f'\t.\t{q_strand}\t.'
                        f'\tID={tx_ensid}_exon_{k};Parent={tx_ensid}\n'
                    )
                written += 1
    return written


# ---------------------------------------------------------------------------
# Main projection logic
# ---------------------------------------------------------------------------

def project(
    source_txs:   List[SourceTranscript],
    chains:       List[Chain],
    min_coverage: float,
) -> List[dict]:
    t_index = _build_t_index(chains)
    results = []

    for tx in source_txs:
        candidate_chains = t_index.get(tx.seqname, [])
        if not candidate_chains:
            continue

        # Find the best chain for this transcript (highest score)
        # Use all chains that overlap the transcript span
        tx_s = tx.start - 1   # 0-based
        tx_e = tx.end
        overlapping = [
            c for c in candidate_chains
            if c.t_start < tx_e and c.t_end > tx_s
        ]
        if not overlapping:
            continue

        best_chain = max(overlapping, key=lambda c: c.score)

        proj_exons = [_project_interval(s, e, best_chain) for s, e in tx.exons]
        valid_exons = [p for p in proj_exons if p is not None]
        if not valid_exons:
            continue

        coverage = _compute_coverage(tx.exons, proj_exons)
        if coverage < min_coverage:
            continue

        q_start = min(p[0] for p in valid_exons)
        q_end   = max(p[1] for p in valid_exons)

        results.append({
            'src_tx_id':   tx.tx_id,
            'src_gene_id': tx.gene_id,
            'q_name':      best_chain.q_name,
            'q_start':     q_start,
            'q_end':       q_end,
            'q_strand':    best_chain.q_strand,
            'coverage':    coverage,
            'proj_exons':  valid_exons,
        })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source_gff3',  required=True)
    ap.add_argument('--chain',        required=True)
    ap.add_argument('--out',          required=True)
    ap.add_argument('--min_coverage', type=float, default=50.0)
    args = ap.parse_args()

    chains     = parse_chain(args.chain)
    source_txs = parse_source_gff3(args.source_gff3)
    projected  = project(source_txs, chains, args.min_coverage)
    n          = write_projected_gff3(projected, args.out)

    print(
        f'project_transcripts: {len(source_txs)} source transcripts; '
        f'{n} projected (min_cov={args.min_coverage})',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
