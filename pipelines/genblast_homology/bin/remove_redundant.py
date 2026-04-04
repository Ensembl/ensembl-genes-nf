#!/usr/bin/env python3
"""
remove_redundant.py — Biotype-priority redundancy removal for GenBlast models.

Reads a merged GFF3, groups overlapping genes on the same chromosome+strand
into clusters, and within each cluster retains only genes of the highest
biotype priority. Within the same tier, all non-overlapping genes are kept.

Priority order (highest first): genblast_1 > genblast_2 > ... > genblast_7

This approximates the Ensembl RemoveRedundantGenes layer-annotation approach:
models of higher-priority biotypes mask lower-priority overlapping models.

Usage:
    remove_redundant.py --gff3 merged.gff3 --out final.gff3
                        [--sample-id SAMPLE]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional

# Priority: lower number = higher priority
BIOTYPE_PRIORITY = {f'genblast_{i}': i for i in range(1, 8)}


def biotype_rank(bt: str) -> int:
    return BIOTYPE_PRIORITY.get(bt, 99)


# ---------------------------------------------------------------------------
# Data structures (reuse same Gene pattern as cluster_igtr.py)
# ---------------------------------------------------------------------------

@dataclass
class Gene:
    seqname: str
    strand: str
    start: int
    end: int
    gene_id: str
    name: str
    biotype: str
    score: str
    lines: list = field(default_factory=list)

    def overlaps(self, other: 'Gene') -> bool:
        return (self.seqname == other.seqname
                and self.strand == other.strand
                and self.start <= other.end
                and self.end >= other.start)


def parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if '=' in part:
            k, _, v = part.partition('=')
            attrs[k.strip()] = v.strip()
    return attrs


def load_genes(path: str) -> list[Gene]:
    genes: list[Gene] = []
    current: Optional[Gene] = None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            attrs = parse_attrs(cols[8])
            if feature == 'gene':
                if current:
                    genes.append(current)
                current = Gene(
                    seqname=cols[0], strand=cols[6],
                    start=int(cols[3]), end=int(cols[4]),
                    gene_id=attrs.get('ID', ''),
                    name=attrs.get('Name', ''),
                    biotype=attrs.get('biotype', 'genblast_7'),
                    score=cols[5],
                    lines=[line],
                )
            elif current:
                current.lines.append(line)
    if current:
        genes.append(current)
    return genes


# ---------------------------------------------------------------------------
# Clustering — same single-linkage sweep as cluster_igtr.py
# ---------------------------------------------------------------------------

def cluster_genes(genes: list[Gene]) -> list[list[Gene]]:
    sorted_genes = sorted(genes, key=lambda g: (g.seqname, g.strand, g.start))
    clusters: list[list[Gene]] = []
    current_cluster: list[Gene] = []
    cluster_end = -1
    cluster_seq = cluster_strand = ''

    for gene in sorted_genes:
        if (gene.seqname != cluster_seq
                or gene.strand != cluster_strand
                or gene.start > cluster_end):
            if current_cluster:
                clusters.append(current_cluster)
            current_cluster = [gene]
            cluster_end = gene.end
            cluster_seq = gene.seqname
            cluster_strand = gene.strand
        else:
            current_cluster.append(gene)
            cluster_end = max(cluster_end, gene.end)

    if current_cluster:
        clusters.append(current_cluster)
    return clusters


# ---------------------------------------------------------------------------
# Biotype-priority selection within each cluster
# ---------------------------------------------------------------------------

def select_from_cluster(cluster: list[Gene]) -> list[Gene]:
    """
    Keep only genes of the highest-priority biotype in this cluster.
    Within the same biotype, keep all (they may not overlap each other).
    """
    best_rank = min(biotype_rank(g.biotype) for g in cluster)
    return [g for g in cluster if biotype_rank(g.biotype) == best_rank]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(clusters: list[list[Gene]], out_path: str) -> tuple[int, int]:
    total_in = sum(len(c) for c in clusters)
    total_out = 0
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for cluster in clusters:
            kept = select_from_cluster(cluster)
            for gene in kept:
                for line in gene.lines:
                    fh.write(line + '\n')
            total_out += len(kept)
    return total_out, total_in - total_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff3',      required=True)
    parser.add_argument('--out',       required=True)
    parser.add_argument('--sample-id', default='sample')
    args = parser.parse_args()

    print(f'[remove_redundant] Loading {args.gff3}', flush=True)
    genes = load_genes(args.gff3)
    print(f'[remove_redundant] {len(genes)} genes loaded', flush=True)

    clusters = cluster_genes(genes)
    print(f'[remove_redundant] {len(clusters)} overlap clusters', flush=True)

    kept, removed = write_output(clusters, args.out)
    print(f'[remove_redundant] Wrote {kept} genes ({removed} removed by priority filter)',
          flush=True)


if __name__ == '__main__':
    main()
