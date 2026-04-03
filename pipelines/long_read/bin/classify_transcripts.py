#!/usr/bin/env python3
"""
classify_transcripts.py — Assign biotypes to collapsed transcript models
based on blastp protein support.

Reads the collapsed GFF3 and a blastp TSV (outfmt 6) and annotates each
transcript with one of:
  isoseq_supported — has a blastp hit at or below the e-value threshold
  isoseq           — no blastp support

The gene biotype is updated to 'isoseq_supported' if ANY transcript has support.

Usage:
    classify_transcripts.py --gff3 in.gff3 --blast hits.tsv --out out.gff3
                            [--evalue 1e-5] [--sample-id SAMPLE]
"""

import argparse
import sys
import re
from pathlib import Path


# ---------------------------------------------------------------------------
# Parse blastp outfmt 6 TSV
# ---------------------------------------------------------------------------

def parse_blast(blast_tsv: str, evalue_cutoff: float) -> set[str]:
    """Return set of query IDs with a hit at or below evalue_cutoff."""
    supported = set()
    with open(blast_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 11:
                continue
            query   = cols[0]
            evalue  = float(cols[10])
            if evalue <= evalue_cutoff:
                supported.add(query)
    return supported


# ---------------------------------------------------------------------------
# Parse + update GFF3
# ---------------------------------------------------------------------------

def classify_gff3(gff3_in: str, gff3_out: str, supported_ids: set[str]) -> dict:
    """
    Read GFF3, update biotype attributes, return counts.
    Transcripts with an ID in supported_ids get biotype=isoseq_supported.
    Genes get biotype=isoseq_supported if any child transcript is supported.
    """
    gene_support: dict[str, bool] = {}   # gene_id -> has_support
    lines_out: list[str] = []
    counts = {'isoseq_supported': 0, 'isoseq': 0, 'genes_supported': 0}

    with open(gff3_in) as fh:
        for line in fh:
            if line.startswith('#'):
                lines_out.append(line)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                lines_out.append(line)
                continue

            feat_type = cols[2]
            attrs     = cols[8]

            if feat_type == 'transcript':
                tx_id   = _get_attr(attrs, 'ID')
                gene_id = _get_attr(attrs, 'Parent')
                if tx_id in supported_ids:
                    attrs = _set_attr(attrs, 'biotype', 'isoseq_supported')
                    gene_support[gene_id] = True
                    counts['isoseq_supported'] += 1
                else:
                    attrs = _set_attr(attrs, 'biotype', 'isoseq')
                    gene_support.setdefault(gene_id, False)
                    counts['isoseq'] += 1
                cols[8] = attrs

            lines_out.append('\t'.join(cols) + '\n')

    # Second pass: update gene biotypes
    final_lines = []
    for line in lines_out:
        if line.startswith('#') or '\t' not in line:
            final_lines.append(line)
            continue
        cols = line.rstrip('\n').split('\t')
        if cols[2] == 'gene':
            gene_id = _get_attr(cols[8], 'ID')
            if gene_support.get(gene_id, False):
                cols[8] = _set_attr(cols[8], 'biotype', 'isoseq_supported')
                counts['genes_supported'] += 1
        final_lines.append('\t'.join(cols) + '\n')

    with open(gff3_out, 'w') as fh:
        fh.writelines(final_lines)

    return counts


def _get_attr(attrs: str, key: str) -> str:
    m = re.search(rf'(?:^|;){key}=([^;]+)', attrs)
    return m.group(1) if m else ''


def _set_attr(attrs: str, key: str, value: str) -> str:
    if re.search(rf'(?:^|;){key}=', attrs):
        return re.sub(rf'((?:^|;)){key}=[^;]*', rf'\g<1>{key}={value}', attrs)
    return attrs + f';{key}={value}'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gff3',      required=True, help='Input collapsed GFF3')
    parser.add_argument('--blast',     required=True, help='blastp outfmt 6 TSV')
    parser.add_argument('--out',       required=True, help='Output classified GFF3')
    parser.add_argument('--evalue',    type=float, default=1e-5,
                        help='E-value cutoff for protein support (default: 1e-5)')
    parser.add_argument('--sample-id', default='sample',
                        help='Sample identifier for logging')
    args = parser.parse_args()

    print(f"[classify_transcripts] Sample: {args.sample_id}", flush=True)
    supported = parse_blast(args.blast, args.evalue)
    print(f"[classify_transcripts] {len(supported)} transcripts with blastp support "
          f"(evalue <= {args.evalue})", flush=True)

    counts = classify_gff3(args.gff3, args.out, supported)
    print(f"[classify_transcripts] isoseq_supported={counts['isoseq_supported']} "
          f"isoseq={counts['isoseq']} genes_with_support={counts['genes_supported']}", flush=True)
    print(f"[classify_transcripts] Written to {args.out}", flush=True)


if __name__ == '__main__':
    main()
