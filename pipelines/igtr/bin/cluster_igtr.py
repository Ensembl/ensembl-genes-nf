#!/usr/bin/env python3
"""
cluster_igtr.py — Cluster overlapping IGTR gene models and select the best.

Reads a merged GFF3 from all GenBlast batches (gene/transcript/exon features),
groups genes that overlap on the same chromosome+strand, and within each
cluster retains the single gene with the highest combined PID+Coverage score.

This emulates the logic of HiveCollapseIGTR: best gene = max(PID + Coverage).

Usage:
    cluster_igtr.py --gff3 merged.gff3 --out clustered.gff3 [--sample-id SAMPLE]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Data structures
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
    # transcript attributes
    tx_id: str = ''
    pid: float = 0.0
    coverage: float = 0.0
    # raw GFF3 lines for this gene block (gene + transcript + exons)
    lines: list = field(default_factory=list)

    @property
    def combined_score(self) -> float:
        return self.pid + self.coverage

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


# ---------------------------------------------------------------------------
# GFF3 parsing — read gene blocks
# ---------------------------------------------------------------------------

def load_genes(gff3_path: str) -> list[Gene]:
    """
    Parse GFF3 into Gene objects.  Every gene feature starts a new block;
    the immediately following transcript line gives us PID and Coverage.
    """
    genes: list[Gene] = []
    current_gene: Optional[Gene] = None

    with open(gff3_path) as fh:
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
                if current_gene:
                    genes.append(current_gene)
                current_gene = Gene(
                    seqname=cols[0],
                    strand=cols[6],
                    start=int(cols[3]),
                    end=int(cols[4]),
                    gene_id=attrs.get('ID', ''),
                    name=attrs.get('Name', ''),
                    biotype=attrs.get('biotype', 'ig_gene'),
                    score=cols[5],
                    lines=[line],
                )

            elif feature == 'transcript' and current_gene:
                try:
                    current_gene.pid      = float(attrs.get('PID', 0))
                    current_gene.coverage = float(attrs.get('Coverage', 0))
                except ValueError:
                    pass
                current_gene.tx_id = attrs.get('ID', '')
                current_gene.lines.append(line)

            elif feature == 'exon' and current_gene:
                current_gene.lines.append(line)

    if current_gene:
        genes.append(current_gene)

    return genes


# ---------------------------------------------------------------------------
# Clustering — single-linkage sweep
# ---------------------------------------------------------------------------

def cluster_genes(genes: list[Gene]) -> list[list[Gene]]:
    """
    Group genes into overlap clusters using a coordinate sweep.
    Sort by (seqname, strand, start); merge into current cluster while
    any gene in the cluster overlaps the next candidate.
    """
    # Sort so we can sweep linearly
    sorted_genes = sorted(genes, key=lambda g: (g.seqname, g.strand, g.start))

    clusters: list[list[Gene]] = []
    current_cluster: list[Gene] = []
    cluster_end = -1
    cluster_seq = ''
    cluster_strand = ''

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
# Output
# ---------------------------------------------------------------------------

def write_clustered(clusters: list[list[Gene]], out_path: str) -> tuple[int, int]:
    """
    For each cluster pick the gene with the highest PID+Coverage and write it.
    Returns (n_clusters, n_removed).
    """
    total_in = sum(len(c) for c in clusters)
    written = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for cluster in clusters:
            best = max(cluster, key=lambda g: g.combined_score)
            for line in best.lines:
                fh.write(line + '\n')
            written += 1

    return written, total_in - written


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff3', required=True, help='Merged GFF3 from all batches')
    parser.add_argument('--out',  required=True, help='Output clustered GFF3')
    parser.add_argument('--sample-id', default='sample', help='Sample ID (informational)')
    args = parser.parse_args()

    print(f'[cluster_igtr] Loading genes from {args.gff3}', flush=True)
    genes = load_genes(args.gff3)
    print(f'[cluster_igtr] Loaded {len(genes)} genes', flush=True)

    clusters = cluster_genes(genes)
    print(f'[cluster_igtr] Found {len(clusters)} overlap clusters', flush=True)

    kept, removed = write_clustered(clusters, args.out)
    print(f'[cluster_igtr] Wrote {kept} genes to {args.out} '
          f'({removed} redundant models removed)', flush=True)


if __name__ == '__main__':
    main()
