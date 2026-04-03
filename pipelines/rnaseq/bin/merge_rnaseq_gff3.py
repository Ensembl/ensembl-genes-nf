#!/usr/bin/env python3
"""
merge_rnaseq_gff3.py — Merge per-sample RNA-seq GFF3 files.

Concatenates GFF3 from multiple samples, re-clusters overlapping transcripts
across samples, and assigns unified gene IDs.  Does NOT deduplicate identical
transcripts — that is handled by downstream gene-building steps.

Biotype is set to 'rnaseq_merged' on all gene features.

Usage:
    merge_rnaseq_gff3.py \\
        --gff3  sample1.gff3 sample2.gff3 ... \\
        --out   merged_rnaseq.gff3
"""

import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List


# ---------------------------------------------------------------------------
# Parsing (re-use pattern from other pipelines)
# ---------------------------------------------------------------------------

@dataclass
class TxRecord:
    seqname: str
    start:   int
    end:     int
    strand:  str
    lines:   List[str] = field(default_factory=list)  # all GFF3 lines for this tx block


def _open(path: str):
    import gzip
    return gzip.open if path.endswith('.gz') else open


def parse_gff3_transcripts(path: str) -> List[TxRecord]:
    """Collect transcript-level blocks (gene+transcript+exon lines)."""
    records: List[TxRecord] = []
    current: TxRecord | None = None

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
            if feature == 'transcript':
                if current is not None:
                    records.append(current)
                current = TxRecord(
                    seqname = cols[0],
                    start   = int(cols[3]),
                    end     = int(cols[4]),
                    strand  = cols[6],
                    lines   = [],
                )
                # Omit the gene line — we'll regenerate it after clustering
                current.lines.append(line)
            elif feature == 'exon' and current is not None:
                current.lines.append(line)
            # skip 'gene' lines — regenerated after clustering

    if current is not None:
        records.append(current)
    return records


# ---------------------------------------------------------------------------
# Clustering (same single-linkage as select_best_targeted)
# ---------------------------------------------------------------------------

def cluster(txs: List[TxRecord]) -> List[List[TxRecord]]:
    if not txs:
        return []
    sorted_txs = sorted(txs, key=lambda t: (t.seqname, t.strand, t.start))
    clusters: List[List[TxRecord]] = []
    cur = [sorted_txs[0]]
    cur_end = sorted_txs[0].end

    for tx in sorted_txs[1:]:
        if (tx.seqname == cur[0].seqname
                and tx.strand == cur[0].strand
                and tx.start <= cur_end):
            cur.append(tx)
            cur_end = max(cur_end, tx.end)
        else:
            clusters.append(cur)
            cur = [tx]
            cur_end = tx.end

    clusters.append(cur)
    return clusters


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def _replace_attr(attrs: str, key: str, value: str) -> str:
    parts = attrs.split(';')
    new_parts = []
    replaced = False
    for p in parts:
        if p.startswith(key + '='):
            new_parts.append(f'{key}={value}')
            replaced = True
        else:
            new_parts.append(p)
    if not replaced:
        new_parts.append(f'{key}={value}')
    return ';'.join(new_parts)


def write_merged_gff3(
    clusters: List[List[TxRecord]],
    out_path: str,
) -> tuple:
    gene_count = 0
    tx_count   = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i, cluster_txs in enumerate(clusters, 1):
            seqname = cluster_txs[0].seqname
            strand  = cluster_txs[0].strand
            g_start = min(tx.start for tx in cluster_txs)
            g_end   = max(tx.end   for tx in cluster_txs)
            gene_id = f'rnaseq_merged_gene_{i:08d}'

            fh.write(
                f'{seqname}\tStringTie2\tgene\t{g_start}\t{g_end}'
                f'\t.\t{strand}\t.'
                f'\tID={gene_id};biotype=rnaseq_merged\n'
            )
            gene_count += 1

            for j, tx in enumerate(sorted(cluster_txs, key=lambda t: t.start), 1):
                tx_id = f'{gene_id}_tx_{j:04d}'
                for k, line in enumerate(tx.lines):
                    cols = line.split('\t')
                    if cols[2] == 'transcript':
                        attrs = _replace_attr(cols[8], 'ID', tx_id)
                        attrs = _replace_attr(attrs, 'Parent', gene_id)
                        attrs = _replace_attr(attrs, 'biotype', 'rnaseq_merged')
                        cols[8] = attrs
                        fh.write('\t'.join(cols) + '\n')
                    elif cols[2] == 'exon':
                        attrs = _replace_attr(cols[8], 'Parent', tx_id)
                        exon_id = f'{tx_id}_exon_{k}'
                        attrs = _replace_attr(attrs, 'ID', exon_id)
                        cols[8] = attrs
                        fh.write('\t'.join(cols) + '\n')
                tx_count += 1

    return gene_count, tx_count


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff3', nargs='+', required=True)
    ap.add_argument('--out',  required=True)
    args = ap.parse_args()

    all_txs: List[TxRecord] = []
    for path in args.gff3:
        all_txs.extend(parse_gff3_transcripts(path))

    clusters      = cluster(all_txs)
    gene_count, tx_count = write_merged_gff3(clusters, args.out)

    print(f'merge_rnaseq_gff3: {len(all_txs)} transcripts from {len(args.gff3)} files; '
          f'{gene_count} merged genes, {tx_count} transcripts written', file=sys.stderr)


if __name__ == '__main__':
    main()
