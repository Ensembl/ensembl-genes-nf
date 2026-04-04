#!/usr/bin/env python3
"""
merge_ab_initio_gff3.py — Merge per-chunk Augustus GFF3 files and renumber IDs.

Each genome chunk produces a separate GFF3 from parse_augustus_gff.py.
This script concatenates them, renumbers all IDs sequentially, and ensures
Parent references are updated to match.

Usage:
    merge_ab_initio_gff3.py \\
        --inputs chunk1.gff3 chunk2.gff3 ... \\
        --out    ab_initio.merged.gff3
"""

import argparse
import re
import sys
from typing import List


def _update_ids(lines: List[str], gene_offset: int) -> List[str]:
    """
    Renumber gene/transcript/exon/CDS IDs in a block of GFF3 lines.
    IDs have the form:  ab_initio_gene_NNNNNNNN[_tx_NNNN[_exon/cds_N]]
    Replace the gene index with gene_offset + original_index.
    """
    # Build a mapping from old gene id → new gene id for this block
    id_map: dict = {}
    result = []

    for line in lines:
        if not line or line.startswith('#'):
            result.append(line)
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            result.append(line)
            continue

        attrs = cols[8]

        # Extract and remap all ab_initio_gene_NNNNNNNN tokens
        def _remap(m):
            old_full = m.group(0)
            gene_num = int(m.group(1))
            new_num  = gene_num + gene_offset
            new_full = old_full.replace(f'gene_{gene_num:08d}', f'gene_{new_num:08d}')
            return new_full

        new_attrs = re.sub(r'ab_initio_gene_(\d{8})', _remap, attrs)
        cols[8] = new_attrs
        result.append('\t'.join(cols))

    return result


def count_genes(lines: List[str]) -> int:
    return sum(1 for l in lines if '\tgene\t' in l and not l.startswith('#'))


def merge(inputs: List[str], out_path: str) -> int:
    gene_offset = 0
    total_genes = 0

    with open(out_path, 'w') as fh:
        fh.write('##gff-version 3\n')
        for path in inputs:
            with open(path) as src:
                lines = [l.rstrip('\n') for l in src if not l.startswith('##')]

            updated = _update_ids(lines, gene_offset)
            chunk_genes = count_genes(updated)
            for line in updated:
                if line:
                    fh.write(line + '\n')
            gene_offset += chunk_genes
            total_genes += chunk_genes

    return total_genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--out',    required=True)
    args = ap.parse_args()

    n = merge(args.inputs, args.out)
    print(f'merge_ab_initio_gff3: {n} genes written to {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
