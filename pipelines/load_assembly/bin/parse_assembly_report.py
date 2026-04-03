#!/usr/bin/env python3
"""
parse_assembly_report.py — Parse NCBI assembly_report.txt and produce:
  1. A seq-region synonyms TSV for coordinate mapping (RefSeq → chr name)
  2. A JSON metadata file with assembly info

NCBI assembly_report.txt column layout (tab-separated):
  0  Sequence-Name          (e.g. 1, X, MT, HSCHR1_CTG1_UNLOCALIZED)
  1  Sequence-Role          (assembled-molecule, unlocalized-scaffold, etc.)
  2  Assigned-Molecule      (e.g. 1, X, MT, na)
  3  Assigned-Molecule-Location/Type  (Chromosome, Mitochondrion, na)
  4  GenBank-Accn           (e.g. CM000663.2)
  5  Relationship           (=, <>)
  6  RefSeq-Accn            (e.g. NC_000001.11)
  7  Assembly-Unit          (Primary Assembly, non-nuclear)
  8  Sequence-Length
  9  UCSC-style-name        (e.g. chr1, chrX, chrM)

Usage:
    parse_assembly_report.py \\
        --report assembly_report.txt \\
        --synonyms_out synonyms.tsv \\
        --meta_out    assembly_meta.json
"""

import argparse
import json
import re
import sys
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Assembly report parsing
# ---------------------------------------------------------------------------

def _open(path: str):
    import gzip
    if path.endswith('.gz'):
        return gzip.open
    return open


def parse_assembly_report(path: str) -> Tuple[List[dict], Dict[str, str]]:
    """
    Parse assembly_report.txt.

    Returns:
      rows     — list of dicts, one per sequence
      metadata — dict of assembly-level key-value pairs from header
    """
    rows = []
    metadata: Dict[str, str] = {}

    opener = _open(path)
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            # Header lines start with #
            if line.startswith('#'):
                # Metadata lines: "# Assembly name:  GRCh38.p14"
                m = re.match(r'^#\s+([^:]+):\s+(.+)$', line)
                if m:
                    key = m.group(1).strip().lower().replace(' ', '_').replace('/', '_')
                    val = m.group(2).strip()
                    metadata[key] = val
                continue

            cols = line.split('\t')
            if len(cols) < 9:
                continue

            row = {
                'seq_name':  cols[0],
                'role':      cols[1],
                'molecule':  cols[2],
                'loc_type':  cols[3],
                'genbank':   cols[4] if cols[4] != 'na' else None,
                'refseq':    cols[6] if cols[6] != 'na' else None,
                'unit':      cols[7],
                'length':    int(cols[8]) if cols[8].isdigit() else 0,
                'ucsc':      cols[9].strip() if len(cols) > 9 and cols[9].strip() not in ('na', '') else None,
            }
            rows.append(row)

    return rows, metadata


# ---------------------------------------------------------------------------
# Synonyms TSV generation
# ---------------------------------------------------------------------------

def write_synonyms(rows: List[dict], out_path: str) -> int:
    """
    Write RefSeq accession → sequence name synonyms TSV.

    Format: refseq_accession<TAB>seq_region_name

    The seq_region_name preferred order: UCSC name (e.g. chr1) → GenBank name → Sequence-Name.
    Rows with no RefSeq accession are skipped.

    Returns number of synonym pairs written.
    """
    written = 0
    with open(out_path, 'w') as fh:
        fh.write('# RefSeq_accession\tseq_region_name\trole\tlength\n')
        for row in rows:
            refseq = row['refseq']
            if not refseq:
                continue
            # Prefer UCSC name, then GenBank, then raw Sequence-Name
            name = row['ucsc'] or row['genbank'] or row['seq_name']
            fh.write(f'{refseq}\t{name}\t{row["role"]}\t{row["length"]}\n')
            written += 1
    return written


def write_metadata(metadata: Dict[str, str], rows: List[dict], out_path: str) -> None:
    """Write assembly metadata JSON."""
    # Summarise sequence counts by role
    role_counts: Dict[str, int] = {}
    total_length = 0
    for row in rows:
        role_counts[row['role']] = role_counts.get(row['role'], 0) + 1
        total_length += row['length']

    doc = {
        'assembly_metadata': metadata,
        'sequence_summary': {
            'total_sequences': len(rows),
            'total_length':    total_length,
            'by_role':         role_counts,
        },
    }
    with open(out_path, 'w') as fh:
        json.dump(doc, fh, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--report',       required=True, help='NCBI assembly_report.txt (.gz ok)')
    ap.add_argument('--synonyms_out', required=True, help='Output synonyms TSV path')
    ap.add_argument('--meta_out',     required=True, help='Output assembly metadata JSON path')
    args = ap.parse_args()

    rows, metadata = parse_assembly_report(args.report)
    n = write_synonyms(rows, args.synonyms_out)
    write_metadata(metadata, rows, args.meta_out)

    print(f'parse_assembly_report: {len(rows)} sequences, {n} RefSeq synonyms written',
          file=sys.stderr)


if __name__ == '__main__':
    main()
