#!/usr/bin/env python3
"""
consolidate_genes.py — Merge gene models from multiple annotation sources.

Implements a layer-annotation strategy: each input GFF3 has an assigned
priority (lower = higher priority).  Overlapping transcripts are clustered
across all sources.  For each locus, transcripts from higher-priority layers
are kept; lower-priority transcripts that overlap a retained model are
suppressed.  Transcripts that do NOT overlap any higher-priority model are
always kept regardless of layer.

Layer priority (lower integer = higher priority):
  0  long_read        IsoSeq / long-read direct evidence
  1  best_targeted    cDNA/protein best-targeted
  2  rnaseq           short-read assembled
  3  projection       cross-genome projection
  4  refseq           RefSeq/external annotation
  5  ab_initio        Augustus ab initio
  6  genblast         genblast homology

The priority order and biotype labels are configurable via the command line.

Usage:
    consolidate_genes.py \\
        --inputs  long_read.gff3:0 rnaseq.gff3:2 ab_initio.gff3:5 \\
        --out     consolidated.gff3

    Each --inputs entry is  <path>:<priority>  (integer priority).
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Transcript:
    seqname: str
    start:   int    # 1-based
    end:     int
    strand:  str
    tx_id:   str
    gene_id: str
    biotype: str
    priority: int
    lines:   List[str] = field(default_factory=list)  # raw GFF3 lines


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _open(path: str):
    return gzip.open if path.endswith('.gz') else open


def _parse_attrs(attr_str: str) -> dict:
    return {m.group(1): m.group(2)
            for m in re.finditer(r'([^=;]+)=([^;]*)', attr_str)}


def parse_gff3(path: str, priority: int) -> List[Transcript]:
    """
    Read transcript + exon lines from a GFF3.
    Gene lines are dropped (will be regenerated during output).
    """
    txs: Dict[str, Transcript] = {}
    exon_lines: Dict[str, List[str]] = defaultdict(list)

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for raw in fh:
            line = raw.rstrip()
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
                    txs[tx_id] = Transcript(
                        seqname  = cols[0],
                        start    = int(cols[3]),
                        end      = int(cols[4]),
                        strand   = cols[6],
                        tx_id    = tx_id,
                        gene_id  = gene_id,
                        biotype  = attrs.get('biotype', 'unknown'),
                        priority = priority,
                        lines    = [line],
                    )
            elif feature in ('exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR'):
                parent = attrs.get('Parent', '')
                if parent:
                    exon_lines[parent].append(line)

    for tx_id, tx in txs.items():
        tx.lines += exon_lines.get(tx_id, [])

    return list(txs.values())


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

def cluster(txs: List[Transcript]) -> List[List[Transcript]]:
    """
    Single-linkage overlap clustering by seqname + strand.
    Two transcripts are linked if they share any genomic overlap.
    """
    if not txs:
        return []

    # Sort by seqname, strand, start
    sorted_txs = sorted(txs, key=lambda t: (t.seqname, t.strand, t.start))

    clusters: List[List[Transcript]] = []
    current_cluster: List[Transcript] = [sorted_txs[0]]
    current_end = sorted_txs[0].end

    for tx in sorted_txs[1:]:
        prev = current_cluster[0]
        if tx.seqname == prev.seqname and tx.strand == prev.strand and tx.start <= current_end:
            current_cluster.append(tx)
            current_end = max(current_end, tx.end)
        else:
            clusters.append(current_cluster)
            current_cluster = [tx]
            current_end = tx.end

    clusters.append(current_cluster)
    return clusters


# ---------------------------------------------------------------------------
# Layer selection within a cluster
# ---------------------------------------------------------------------------

def select_from_cluster(cluster_txs: List[Transcript]) -> List[Transcript]:
    """
    For one genomic locus:
    1. Sort transcripts by priority (ascending = higher priority first).
    2. Keep all transcripts from the highest-priority layer present.
    3. For each lower-priority transcript, suppress it if it overlaps
       any already-retained transcript; otherwise keep it.
    """
    if not cluster_txs:
        return []

    best_priority = min(t.priority for t in cluster_txs)
    retained: List[Transcript] = []

    # Process priority layers in order
    for priority in sorted({t.priority for t in cluster_txs}):
        layer_txs = [t for t in cluster_txs if t.priority == priority]
        for tx in layer_txs:
            if priority == best_priority:
                retained.append(tx)
            else:
                # Keep only if no overlap with any retained transcript
                overlap = any(
                    tx.seqname == r.seqname
                    and tx.strand == r.strand
                    and tx.start <= r.end
                    and tx.end >= r.start
                    for r in retained
                )
                if not overlap:
                    retained.append(tx)

    return retained


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def _replace_attr(attr_str: str, key: str, value: str) -> str:
    """Replace or add a key=value in a GFF3 attribute string."""
    if re.search(rf'\b{key}=', attr_str):
        return re.sub(rf'\b{key}=[^;]*', f'{key}={value}', attr_str)
    return attr_str.rstrip(';') + f';{key}={value}'


def write_consolidated_gff3(
    projections: List[List[Transcript]],  # list of per-locus selected txs
    out_path: str,
) -> Tuple[int, int]:
    """Write consolidated GFF3. Returns (gene_count, tx_count)."""
    # Flatten and group by source gene_id to regenerate gene features
    all_txs = [tx for locus in projections for tx in locus]
    by_gene: Dict[str, List[Transcript]] = defaultdict(list)
    for tx in all_txs:
        by_gene[tx.gene_id].append(tx)

    gene_count = tx_count = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, (orig_gene_id, txs) in enumerate(by_gene.items(), 1):
            seqname  = txs[0].seqname
            strand   = txs[0].strand
            g_start  = min(t.start for t in txs)
            g_end    = max(t.end   for t in txs)
            biotype  = txs[0].biotype
            gene_id  = f'consolidated_gene_{i:08d}'

            fh.write(
                f'{seqname}\tconsolidation\tgene'
                f'\t{g_start}\t{g_end}\t.\t{strand}\t.'
                f'\tID={gene_id};Name={orig_gene_id};biotype={biotype};'
                f'src_gene={orig_gene_id}\n'
            )
            gene_count += 1

            for tx in txs:
                # Rewrite transcript line with new Parent
                tx_lines = list(tx.lines)
                if tx_lines:
                    cols = tx_lines[0].split('\t')
                    if len(cols) >= 9:
                        new_attrs = _replace_attr(cols[8], 'Parent', gene_id)
                        cols[8] = new_attrs
                        tx_lines[0] = '\t'.join(cols)
                for line in tx_lines:
                    fh.write(line + '\n')
                tx_count += 1

    return gene_count, tx_count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--inputs', nargs='+', required=True,
        metavar='PATH:PRIORITY',
        help='GFF3 files with priority: e.g. rnaseq.gff3:2'
    )
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    all_txs: List[Transcript] = []
    for spec in args.inputs:
        if ':' not in spec:
            print(f'ERROR: --inputs must be PATH:PRIORITY, got {spec!r}', file=sys.stderr)
            sys.exit(1)
        path, prio_str = spec.rsplit(':', 1)
        priority = int(prio_str)
        txs = parse_gff3(path, priority)
        all_txs.extend(txs)
        print(f'  Loaded {len(txs)} transcripts from {path} (priority={priority})', file=sys.stderr)

    clusters = cluster(all_txs)
    selected_per_locus = [select_from_cluster(cl) for cl in clusters]

    g_count, tx_count = write_consolidated_gff3(selected_per_locus, args.out)

    total_in  = len(all_txs)
    total_out = sum(len(s) for s in selected_per_locus)
    print(
        f'consolidate_genes: {total_in} input transcripts → '
        f'{total_out} retained in {g_count} genes '
        f'({len(clusters)} loci)',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
