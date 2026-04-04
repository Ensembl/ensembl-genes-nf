#!/usr/bin/env python3
"""
select_best_targeted.py — Merge cDNA and protein exonerate GFF3 files and
select the best non-redundant set of gene models.

Strategy (mirrors HiveBestTargetted logic):
  1. Cluster overlapping transcripts on the same strand/seqname.
  2. For single-analysis clusters: keep all models.
  3. For multi-analysis clusters: prefer cdna_alignment over protein_alignment.
     Within the same biotype, keep the model with the highest (coverage + pid).
  4. Emit merged GFF3 with final biotype 'best_targeted'.

Usage:
    select_best_targeted.py \\
        --cdna_gff3    merged_cdna.gff3 \\
        --protein_gff3 merged_protein.gff3 \\
        --out          best_targeted.gff3 \\
        --min_coverage 50 \\
        --min_pid      50
"""

import argparse
import sys
from dataclasses import dataclass, field
from typing import List, Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Gene:
    seqname:  str
    start:    int
    end:      int
    strand:   str
    score:    float
    coverage: float
    pid:      float
    query_id: str
    biotype:  str
    raw_lines: List[str] = field(default_factory=list)  # original GFF3 lines for this gene block

    @property
    def quality(self) -> float:
        return self.coverage + self.pid

    @property
    def analysis_priority(self) -> int:
        """Lower is better: cdna first."""
        return 0 if 'cdna' in self.biotype else 1


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def _attr(attrs_str: str, key: str) -> str:
    for part in attrs_str.split(';'):
        if part.startswith(key + '='):
            return part[len(key) + 1:]
    return ''


def parse_gff3(path: str) -> List[Gene]:
    if not path:
        return []
    genes: List[Gene] = []
    current: Optional[Gene] = None

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
                if current:
                    genes.append(current)
                attrs = cols[8]
                current = Gene(
                    seqname  = cols[0],
                    start    = int(cols[3]),
                    end      = int(cols[4]),
                    strand   = cols[6],
                    score    = float(cols[5]) if cols[5] != '.' else 0.0,
                    coverage = float(_attr(attrs, 'coverage') or '0'),
                    pid      = float(_attr(attrs, 'pid') or '0'),
                    query_id = _attr(attrs, 'Name') or _attr(attrs, 'ID'),
                    biotype  = _attr(attrs, 'biotype') or 'best_targeted',
                    raw_lines = [line],
                )
            elif current is not None:
                current.raw_lines.append(line)

    if current:
        genes.append(current)
    return genes


def _open(path: str):
    import gzip
    if path.endswith('.gz'):
        return gzip.open
    return open


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

def _overlaps(a: Gene, b: Gene) -> bool:
    return (a.seqname == b.seqname
            and a.strand == b.strand
            and a.start <= b.end
            and b.start <= a.end)


def cluster(genes: List[Gene]) -> List[List[Gene]]:
    """Simple single-linkage overlap clustering."""
    if not genes:
        return []
    sorted_genes = sorted(genes, key=lambda g: (g.seqname, g.strand, g.start))
    clusters: List[List[Gene]] = []
    current_cluster: List[Gene] = [sorted_genes[0]]
    cluster_end = sorted_genes[0].end

    for gene in sorted_genes[1:]:
        if (gene.seqname == current_cluster[0].seqname
                and gene.strand == current_cluster[0].strand
                and gene.start <= cluster_end):
            current_cluster.append(gene)
            cluster_end = max(cluster_end, gene.end)
        else:
            clusters.append(current_cluster)
            current_cluster = [gene]
            cluster_end = gene.end

    clusters.append(current_cluster)
    return clusters


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def select_from_cluster(cluster: List[Gene]) -> List[Gene]:
    """
    For a cluster of overlapping genes:
    - If only one analysis type present: keep all.
    - If mixed: keep only the best-priority analysis type,
      then within that type keep only the single highest-quality model.
    """
    biotypes = {g.biotype for g in cluster}
    analyses = {g.analysis_priority for g in cluster}

    if len(analyses) == 1:
        # Single analysis — keep all
        return list(cluster)

    # Mixed — keep best analysis priority only
    best_priority = min(analyses)
    candidates = [g for g in cluster if g.analysis_priority == best_priority]
    # Among candidates, keep single best by quality
    best = max(candidates, key=lambda g: g.quality)
    return [best]


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_best_gff3(genes: List[Gene], out_path: str) -> int:
    """Rewrite selected genes as GFF3 with biotype=best_targeted."""
    written = 0
    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for gene in genes:
            # Re-emit lines, overriding biotype in gene feature
            for line in gene.raw_lines:
                cols = line.split('\t')
                if len(cols) >= 9 and cols[2] == 'gene':
                    # Replace biotype attribute
                    attrs = cols[8]
                    attrs = _replace_attr(attrs, 'biotype', 'best_targeted')
                    cols[8] = attrs
                    fh.write('\t'.join(cols) + '\n')
                else:
                    fh.write(line + '\n')
            written += 1
    return written


def _replace_attr(attrs: str, key: str, value: str) -> str:
    parts = attrs.split(';')
    replaced = False
    new_parts = []
    for p in parts:
        if p.startswith(key + '='):
            new_parts.append(f'{key}={value}')
            replaced = True
        else:
            new_parts.append(p)
    if not replaced:
        new_parts.append(f'{key}={value}')
    return ';'.join(new_parts)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cdna_gff3',    default=None)
    ap.add_argument('--protein_gff3', default=None)
    ap.add_argument('--out',          required=True)
    ap.add_argument('--min_coverage', type=float, default=50.0)
    ap.add_argument('--min_pid',      type=float, default=50.0)
    args = ap.parse_args()

    cdna_genes    = parse_gff3(args.cdna_gff3)    if args.cdna_gff3    else []
    protein_genes = parse_gff3(args.protein_gff3) if args.protein_gff3 else []

    all_genes = cdna_genes + protein_genes

    # Apply final thresholds
    all_genes = [g for g in all_genes
                 if g.coverage >= args.min_coverage and g.pid >= args.min_pid]

    clusters   = cluster(all_genes)
    selected   = []
    for cl in clusters:
        selected.extend(select_from_cluster(cl))

    n = write_best_gff3(selected, args.out)
    print(f'select_best_targeted: {len(cdna_genes)} cdna + {len(protein_genes)} protein genes; '
          f'{n} selected after clustering', file=sys.stderr)


if __name__ == '__main__':
    main()
