#!/usr/bin/env python3
"""
add_utrs.py — Add UTR regions to consolidated coding gene models.

Algorithm
---------
For each transcript in the consolidated (acceptor) GFF3:

  1. Extract its CDS exon coordinates sorted by position.

  2. Search each donor GFF3 file in turn (first file = highest priority).
     Within each file, search all donor transcripts on the same chromosome
     and strand for one whose exon structure CONTAINS the acceptor CDS:

       • All *internal* CDS splice sites (right-end of every CDS exon except
         the last, and left-end of every CDS exon except the first) must appear
         as exact exon boundaries in the donor transcript.
       • The leftmost boundary of the first CDS exon may differ (UTR extends
         further left in the donor).
       • The rightmost boundary of the last CDS exon may differ (UTR extends
         further right in the donor).

     Single-exon CDS (no internal junctions) requires only that the donor
     exon fully spans the CDS.

  3. Use the first matching donor found (priority = file order, then file order
     within the file).

  4. Rebuild the acceptor transcript's exon list:
       • Keep the CDS exons exactly as they are.
       • Prepend any donor exon(s) that lie 5' of the first CDS exon (trimmed
         to meet the first CDS exon start).
       • Append any donor exon(s) that lie 3' of the last CDS exon (trimmed to
         start at the last CDS exon end).

  5. Clip extended UTR exons so the total UTR does not exceed max_5prime_utr
     or max_3prime_utr bp.

  6. Discard any UTR exon whose final size is below min_utr_exon_size bp.

  7. Write the modified acceptors back as GFF3.  Transcripts that received UTR
     from a donor gain the attribute  utr_source=<donor_file_basename>.

GFF3 notes
----------
  • 1-based, fully-closed coordinates throughout.
  • Strand is '+' or '-'.
  • On the '-' strand "5' UTR" is genomically to the *right* of the CDS and
    "3' UTR" is to the *left*, but the logic below is coordinate-based and
    works correctly for both strands without special-casing.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CdsExon:
    """A single CDS interval (1-based, closed)."""
    start: int
    end: int
    phase: str  # GFF3 phase column (0, 1, 2, or '.')


@dataclass
class Exon:
    """A single exon interval (1-based, closed)."""
    start: int
    end: int


@dataclass
class Transcript:
    """
    A single transcript with its child exons and CDS segments.

    Attributes
    ----------
    seqname : chromosome / sequence name
    start   : transcript start (1-based, equals min exon start)
    end     : transcript end   (1-based, equals max exon end)
    strand  : '+' or '-'
    tx_id   : value of the ID= attribute
    gene_id : value of the Parent= attribute (parent gene)
    source  : GFF3 source column (column 2)
    score   : GFF3 score column  (column 6)
    attrs   : raw attribute string from the mRNA/transcript line
    exons   : list of Exon, sorted by start
    cds     : list of CdsExon, sorted by start
    raw_other_lines : non-exon, non-CDS child feature lines (e.g. UTR lines
                      from the original file — these are dropped on output so
                      the freshly built exon set is the source of truth)
    """
    seqname: str
    start: int
    end: int
    strand: str
    tx_id: str
    gene_id: str
    source: str
    score: str
    attrs: str
    exons: List[Exon] = field(default_factory=list)
    cds: List[CdsExon] = field(default_factory=list)
    raw_other_lines: List[str] = field(default_factory=list)


@dataclass
class Gene:
    """A gene feature with its child transcripts."""
    seqname: str
    start: int
    end: int
    strand: str
    gene_id: str
    source: str
    score: str
    attrs: str
    transcripts: List[Transcript] = field(default_factory=list)


# ---------------------------------------------------------------------------
# GFF3 parsing helpers
# ---------------------------------------------------------------------------

def _open_file(path: str):
    """Return an open file handle for plain or gzip-compressed files."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    """
    Parse a GFF3 attribute column into a dict.

    Handles the most common GFF3 attribute syntax:
      key=value;key2=value2
    Values with commas (multi-value) are returned as-is (first value would
    require splitting on ',').
    """
    result: Dict[str, str] = {}
    for part in attr_str.strip().split(';'):
        part = part.strip()
        if '=' in part:
            k, _, v = part.partition('=')
            result[k.strip()] = v.strip()
    return result


def _set_attr(attr_str: str, key: str, value: str) -> str:
    """
    Set or update a single key=value pair in a GFF3 attribute string.
    If the key already exists it is replaced; otherwise it is appended.
    """
    # Remove existing key=value if present
    cleaned = re.sub(rf'(^|;){re.escape(key)}=[^;]*', '', attr_str)
    cleaned = cleaned.strip(';')
    return f'{cleaned};{key}={value}' if cleaned else f'{key}={value}'


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def parse_gff3(path: str) -> Dict[str, Gene]:
    """
    Parse a GFF3 file into a dict of  gene_id → Gene.

    Recognised feature types
    ------------------------
    gene        → creates a Gene record
    mRNA / transcript → creates a Transcript record, attached to parent Gene
    exon        → added to parent Transcript.exons
    CDS         → added to parent Transcript.cds
    Other child features (five_prime_UTR, three_prime_UTR, …) are stored in
    Transcript.raw_other_lines but are *not* emitted on output (the rebuilt
    exon structure replaces them).

    Transcripts without a parent gene (orphan) are wrapped in a synthetic
    Gene to avoid data loss.

    Returns
    -------
    dict mapping gene_id → Gene (insertion order preserved)
    """
    genes: Dict[str, Gene] = {}
    # tx_id → Transcript, populated while reading; genes may appear after txs
    txs: Dict[str, Transcript] = {}
    # deferred: tx_id → list of (feature_type, cols) for exon/CDS lines that
    # arrived before their parent transcript was seen (shouldn't normally
    # happen in valid GFF3, but be defensive)
    deferred_children: Dict[str, List[Tuple[str, List[str]]]] = defaultdict(list)

    with _open_file(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue

            seqname, source, feature = cols[0], cols[1], cols[2]
            start_s, end_s = cols[3], cols[4]
            score, strand, phase = cols[5], cols[6], cols[7]
            attr_str = cols[8]
            attrs = _parse_attrs(attr_str)

            if feature == 'gene':
                gene_id = attrs.get('ID', f'gene_{seqname}_{start_s}_{end_s}')
                genes[gene_id] = Gene(
                    seqname=seqname,
                    start=int(start_s),
                    end=int(end_s),
                    strand=strand,
                    gene_id=gene_id,
                    source=source,
                    score=score,
                    attrs=attr_str,
                )

            elif feature in ('mRNA', 'transcript'):
                tx_id = attrs.get('ID', '')
                gene_id = attrs.get('Parent', '')
                if not tx_id:
                    continue  # malformed — skip
                tx = Transcript(
                    seqname=seqname,
                    start=int(start_s),
                    end=int(end_s),
                    strand=strand,
                    tx_id=tx_id,
                    gene_id=gene_id,
                    source=source,
                    score=score,
                    attrs=attr_str,
                )
                txs[tx_id] = tx
                # Attach deferred children accumulated before this line
                for (ftype, child_cols) in deferred_children.pop(tx_id, []):
                    _attach_child(tx, ftype, child_cols)

            elif feature == 'exon':
                parent_id = attrs.get('Parent', '')
                if parent_id in txs:
                    txs[parent_id].exons.append(
                        Exon(start=int(start_s), end=int(end_s))
                    )
                else:
                    deferred_children[parent_id].append((feature, cols))

            elif feature == 'CDS':
                parent_id = attrs.get('Parent', '')
                if parent_id in txs:
                    txs[parent_id].cds.append(
                        CdsExon(start=int(start_s), end=int(end_s), phase=phase)
                    )
                else:
                    deferred_children[parent_id].append((feature, cols))

            else:
                # Other child features (UTR, start_codon, stop_codon, …)
                parent_id = attrs.get('Parent', '')
                if parent_id and parent_id in txs:
                    txs[parent_id].raw_other_lines.append(line)

    # Handle any remaining deferred children (parent gene appears after child)
    for tx_id, children in deferred_children.items():
        if tx_id in txs:
            for ftype, child_cols in children:
                _attach_child(txs[tx_id], ftype, child_cols)

    # Sort exons / CDS within each transcript
    for tx in txs.values():
        tx.exons.sort(key=lambda e: e.start)
        tx.cds.sort(key=lambda c: c.start)

    # Attach transcripts to their parent genes
    for tx in txs.values():
        if tx.gene_id in genes:
            genes[tx.gene_id].transcripts.append(tx)
        else:
            # Orphan transcript — create a synthetic gene wrapper
            synthetic_gene_id = f'gene__{tx.tx_id}'
            if synthetic_gene_id not in genes:
                genes[synthetic_gene_id] = Gene(
                    seqname=tx.seqname,
                    start=tx.start,
                    end=tx.end,
                    strand=tx.strand,
                    gene_id=synthetic_gene_id,
                    source=tx.source,
                    score='.',
                    attrs=f'ID={synthetic_gene_id}',
                )
            genes[synthetic_gene_id].transcripts.append(tx)

    return genes


def _attach_child(tx: Transcript, feature: str, cols: List[str]) -> None:
    """Attach a deferred exon/CDS/other line to a transcript."""
    start, end, phase = int(cols[3]), int(cols[4]), cols[7]
    if feature == 'exon':
        tx.exons.append(Exon(start=start, end=end))
    elif feature == 'CDS':
        tx.cds.append(CdsExon(start=start, end=end, phase=phase))
    else:
        tx.raw_other_lines.append('\t'.join(cols))


# ---------------------------------------------------------------------------
# UTR matching logic
# ---------------------------------------------------------------------------

def cds_matches_donor(
    cds_exons: List[CdsExon],
    donor_exons: List[Exon],
) -> bool:
    """
    Return True if *donor_exons* fully contains all *cds_exons*.

    Rules
    -----
    1. The donor must have at least one exon.
    2. For multi-exon CDS:
       • Every internal splice site of the CDS must be present as a boundary
         in the donor exon list.
         Internal CDS right-ends  (all except the last CDS exon's right end)
         must appear as a donor exon's right end (end coordinate).
         Internal CDS left-ends   (all except the first CDS exon's left end)
         must appear as a donor exon's left end (start coordinate).
       • The leftmost CDS exon start (5'-edge) may be >= the donor exon's
         start (UTR extends leftward in the donor).
       • The rightmost CDS exon end (3'-edge) may be <= the donor exon's end
         (UTR extends rightward in the donor).
    3. For single-exon CDS:
       • There must be exactly one donor exon that spans the entire CDS
         (donor start ≤ CDS start AND donor end ≥ CDS end).

    Parameters
    ----------
    cds_exons    : sorted list of CdsExon (by start, ascending)
    donor_exons  : sorted list of Exon (by start, ascending)

    Returns
    -------
    bool
    """
    if not cds_exons or not donor_exons:
        return False

    if len(cds_exons) == 1:
        # Single-exon CDS — the donor must contain a single exon spanning it
        cds = cds_exons[0]
        for dex in donor_exons:
            if dex.start <= cds.start and dex.end >= cds.end:
                return True
        return False

    # Build fast lookup sets from donor exon boundaries
    donor_starts = {e.start for e in donor_exons}
    donor_ends   = {e.end   for e in donor_exons}

    # Check all internal right-ends of CDS (all but the last exon)
    for cds_ex in cds_exons[:-1]:
        if cds_ex.end not in donor_ends:
            return False

    # Check all internal left-ends of CDS (all but the first exon)
    for cds_ex in cds_exons[1:]:
        if cds_ex.start not in donor_starts:
            return False

    # Verify that the first CDS exon is covered by some donor exon
    first_cds = cds_exons[0]
    if not any(dex.start <= first_cds.start and dex.end >= first_cds.end
               for dex in donor_exons):
        return False

    # Verify that the last CDS exon is covered by some donor exon
    last_cds = cds_exons[-1]
    if not any(dex.start <= last_cds.start and dex.end >= last_cds.end
               for dex in donor_exons):
        return False

    return True


def find_matching_donor(
    acceptor_tx: Transcript,
    donor_genes_by_chrom_strand: Dict[Tuple[str, str], List[Transcript]],
) -> Optional[Transcript]:
    """
    Return the first donor transcript that matches *acceptor_tx*'s CDS.

    Parameters
    ----------
    acceptor_tx
        The transcript whose CDS we want to decorate with UTRs.
    donor_genes_by_chrom_strand
        Pre-built index: (chrom, strand) → list of donor Transcript objects,
        already sorted by (start, tx_id) for deterministic output.

    Returns
    -------
    Matching Transcript or None.
    """
    if not acceptor_tx.cds:
        return None  # No CDS — nothing to match against

    key = (acceptor_tx.seqname, acceptor_tx.strand)
    candidates = donor_genes_by_chrom_strand.get(key, [])

    for donor_tx in candidates:
        if cds_matches_donor(acceptor_tx.cds, donor_tx.exons):
            return donor_tx

    return None


def _build_donor_index(
    donor_genes: Dict[str, Gene],
) -> Dict[Tuple[str, str], List[Transcript]]:
    """
    Index donor transcripts by (chrom, strand) for fast lookup.

    Returns
    -------
    dict : (chrom, strand) → sorted list of Transcript
    """
    index: Dict[Tuple[str, str], List[Transcript]] = defaultdict(list)
    for gene in donor_genes.values():
        for tx in gene.transcripts:
            key = (tx.seqname, tx.strand)
            index[key].append(tx)
    # Sort for determinism
    for key in index:
        index[key].sort(key=lambda t: (t.start, t.tx_id))
    return index


# ---------------------------------------------------------------------------
# UTR grafting
# ---------------------------------------------------------------------------

def add_utr_to_transcript(
    acceptor_tx: Transcript,
    donor_tx: Transcript,
    max_5prime: int,
    max_3prime: int,
    min_exon: int,
) -> List[Exon]:
    """
    Build a new exon list for *acceptor_tx* by grafting UTR regions from
    *donor_tx* onto the acceptor's CDS exons.

    The CDS exon positions are preserved exactly; only the flanking UTR
    exons are derived from the donor.

    Parameters
    ----------
    acceptor_tx : Transcript with .cds populated and .exons as the current
                  (possibly UTR-less) exon list
    donor_tx    : Transcript whose .exons contain UTR regions
    max_5prime  : maximum 5' UTR length in bp (total, across all 5' UTR exons)
    max_3prime  : maximum 3' UTR length in bp
    min_exon    : minimum UTR exon size in bp to keep

    Returns
    -------
    New list of Exon sorted by start position.
    """
    cds_exons = sorted(acceptor_tx.cds, key=lambda c: c.start)
    cds_start = cds_exons[0].start   # genomic start of CDS
    cds_end   = cds_exons[-1].end    # genomic end   of CDS

    donor_exons = sorted(donor_tx.exons, key=lambda e: e.start)

    # ------------------------------------------------------------------ #
    #  Collect 5' UTR exons (genomically left of cds_start)               #
    # ------------------------------------------------------------------ #
    five_prime_exons: List[Exon] = []
    for dex in donor_exons:
        if dex.end < cds_start:
            # Entirely in the UTR
            five_prime_exons.append(Exon(start=dex.start, end=dex.end))
        elif dex.start < cds_start <= dex.end:
            # Exon straddles the CDS start — take only the UTR portion
            five_prime_exons.append(Exon(start=dex.start, end=cds_start - 1))
        # Exons at/after cds_start are handled below

    # ------------------------------------------------------------------ #
    #  Collect 3' UTR exons (genomically right of cds_end)                #
    # ------------------------------------------------------------------ #
    three_prime_exons: List[Exon] = []
    for dex in donor_exons:
        if dex.start > cds_end:
            # Entirely in the UTR
            three_prime_exons.append(Exon(start=dex.start, end=dex.end))
        elif dex.start <= cds_end < dex.end:
            # Exon straddles the CDS end — take only the UTR portion
            three_prime_exons.append(Exon(start=cds_end + 1, end=dex.end))

    # ------------------------------------------------------------------ #
    #  Apply max-length limits                                            #
    # ------------------------------------------------------------------ #
    five_prime_exons  = _clip_utr(five_prime_exons,  max_5prime,  from_right=True)
    three_prime_exons = _clip_utr(three_prime_exons, max_3prime,  from_right=False)

    # ------------------------------------------------------------------ #
    #  Filter exons smaller than min_exon                                 #
    # ------------------------------------------------------------------ #
    five_prime_exons  = [e for e in five_prime_exons  if (e.end - e.start + 1) >= min_exon]
    three_prime_exons = [e for e in three_prime_exons if (e.end - e.start + 1) >= min_exon]

    # ------------------------------------------------------------------ #
    #  Assemble final exon list: 5'UTR + CDS exons + 3'UTR                #
    # ------------------------------------------------------------------ #
    cds_as_exons = [Exon(start=c.start, end=c.end) for c in cds_exons]
    all_exons = five_prime_exons + cds_as_exons + three_prime_exons
    all_exons.sort(key=lambda e: e.start)
    return all_exons


def _clip_utr(
    utr_exons: List[Exon],
    max_len: int,
    from_right: bool,
) -> List[Exon]:
    """
    Trim a list of UTR exons so their total span does not exceed *max_len* bp.

    Parameters
    ----------
    utr_exons  : list of Exon sorted by start (ascending)
    max_len    : maximum total UTR length in bp
    from_right : if True, trimming removes from the left side of the list
                 (i.e. we keep exons closest to the CDS = right-most for 5' UTR);
                 if False, trimming removes from the right side of the list
                 (i.e. we keep exons closest to the CDS = left-most for 3' UTR).

    Returns
    -------
    Clipped list of Exon (still sorted by start).
    """
    if not utr_exons:
        return []

    total = sum(e.end - e.start + 1 for e in utr_exons)
    if total <= max_len:
        return utr_exons

    clipped: List[Exon] = []
    remaining = max_len

    # Iterate from the CDS-proximal side
    ordered = list(reversed(utr_exons)) if from_right else list(utr_exons)
    for exon in ordered:
        size = exon.end - exon.start + 1
        if size <= remaining:
            clipped.append(exon)
            remaining -= size
        else:
            # Partial exon
            if from_right:
                clipped.append(Exon(start=exon.end - remaining + 1, end=exon.end))
            else:
                clipped.append(Exon(start=exon.start, end=exon.start + remaining - 1))
            remaining = 0
        if remaining <= 0:
            break

    clipped.sort(key=lambda e: e.start)
    return clipped


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def _gff_row(
    seqname: str,
    source: str,
    feature: str,
    start: int,
    end: int,
    score: str,
    strand: str,
    phase: str,
    attrs: str,
) -> str:
    """Format a single GFF3 data row (no trailing newline)."""
    return '\t'.join([
        seqname, source, feature,
        str(start), str(end),
        score, strand, phase,
        attrs,
    ])


def write_gff3(
    genes: Dict[str, Gene],
    output_path: str,
) -> Tuple[int, int]:
    """
    Write all genes (with their updated transcripts) to a GFF3 file.

    Each gene line is regenerated from the Gene object.  Each transcript
    line is regenerated with any updated attributes (e.g. utr_source).
    Each exon line is regenerated from the Transcript.exons list.
    CDS lines are preserved from Transcript.cds.

    Parameters
    ----------
    genes       : dict gene_id → Gene (with updated Transcript.exons)
    output_path : path to write GFF3

    Returns
    -------
    (gene_count, transcript_count)
    """
    gene_count = tx_count = 0

    with open(output_path, 'w') as fh:
        fh.write('##gff-version 3\n')

        for gene in genes.values():
            if not gene.transcripts:
                continue

            # Recalculate gene boundaries from its transcripts
            g_start = min(tx.start for tx in gene.transcripts)
            g_end   = max(tx.end   for tx in gene.transcripts)
            # Adjust for any newly added UTR exons
            for tx in gene.transcripts:
                if tx.exons:
                    g_start = min(g_start, tx.exons[0].start)
                    g_end   = max(g_end,   tx.exons[-1].end)

            fh.write(_gff_row(
                seqname=gene.seqname,
                source=gene.source,
                feature='gene',
                start=g_start,
                end=g_end,
                score=gene.score,
                strand=gene.strand,
                phase='.',
                attrs=gene.attrs,
            ) + '\n')
            gene_count += 1

            for tx in gene.transcripts:
                # Transcript coordinates updated to span all exons
                tx_start = tx.exons[0].start if tx.exons else tx.start
                tx_end   = tx.exons[-1].end  if tx.exons else tx.end

                fh.write(_gff_row(
                    seqname=tx.seqname,
                    source=tx.source,
                    feature='mRNA',
                    start=tx_start,
                    end=tx_end,
                    score=tx.score,
                    strand=tx.strand,
                    phase='.',
                    attrs=tx.attrs,
                ) + '\n')
                tx_count += 1

                # Exon lines
                for exon_num, exon in enumerate(tx.exons, 1):
                    tx_base_id = _parse_attrs(tx.attrs).get('ID', tx.tx_id)
                    fh.write(_gff_row(
                        seqname=tx.seqname,
                        source=tx.source,
                        feature='exon',
                        start=exon.start,
                        end=exon.end,
                        score='.',
                        strand=tx.strand,
                        phase='.',
                        attrs=f'Parent={tx_base_id}',
                    ) + '\n')

                # CDS lines
                for cds_ex in sorted(tx.cds, key=lambda c: c.start):
                    tx_base_id = _parse_attrs(tx.attrs).get('ID', tx.tx_id)
                    fh.write(_gff_row(
                        seqname=tx.seqname,
                        source=tx.source,
                        feature='CDS',
                        start=cds_ex.start,
                        end=cds_ex.end,
                        score='.',
                        strand=tx.strand,
                        phase=cds_ex.phase,
                        attrs=f'Parent={tx_base_id}',
                    ) + '\n')

    return gene_count, tx_count


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description='Add UTR regions to consolidated coding gene models.'
    )
    ap.add_argument('--consolidated', required=True,
                    help='GFF3 of consolidated coding gene models (acceptors)')
    ap.add_argument('--donors', required=True, nargs='+',
                    help='Donor GFF3 files in priority order (first = highest priority)')
    ap.add_argument('--out', required=True,
                    help='Output GFF3 path')
    ap.add_argument('--max-5prime', type=int, default=5000,
                    help='Maximum 5\' UTR length in bp (default: 5000)')
    ap.add_argument('--max-3prime', type=int, default=10000,
                    help='Maximum 3\' UTR length in bp (default: 10000)')
    ap.add_argument('--min-utr-exon', type=int, default=30,
                    help='Minimum UTR exon size in bp to retain (default: 30)')
    args = ap.parse_args()

    # ------------------------------------------------------------------ #
    #  Load acceptor models                                               #
    # ------------------------------------------------------------------ #
    print(f'Loading acceptor models from: {args.consolidated}', file=sys.stderr)
    acceptor_genes = parse_gff3(args.consolidated)
    n_acceptors = sum(len(g.transcripts) for g in acceptor_genes.values())
    print(f'  {len(acceptor_genes)} genes / {n_acceptors} transcripts loaded',
          file=sys.stderr)

    # ------------------------------------------------------------------ #
    #  Load donor models (in priority order) and build indexes           #
    # ------------------------------------------------------------------ #
    donor_indexes: List[Tuple[str, Dict[Tuple[str, str], List[Transcript]]]] = []
    for donor_path in args.donors:
        print(f'Loading donor models from: {donor_path}', file=sys.stderr)
        donor_genes = parse_gff3(donor_path)
        n_donors = sum(len(g.transcripts) for g in donor_genes.values())
        print(f'  {len(donor_genes)} genes / {n_donors} transcripts loaded',
              file=sys.stderr)
        donor_basename = os.path.basename(donor_path)
        index = _build_donor_index(donor_genes)
        donor_indexes.append((donor_basename, index))

    # ------------------------------------------------------------------ #
    #  Match acceptors to donors and graft UTRs                          #
    # ------------------------------------------------------------------ #
    n_matched = 0
    n_skipped_no_cds = 0

    for gene in acceptor_genes.values():
        for tx in gene.transcripts:
            if not tx.cds:
                n_skipped_no_cds += 1
                continue

            matched_donor_tx: Optional[Transcript] = None
            matched_source: str = ''

            # Try each donor file in priority order
            for donor_basename, donor_index in donor_indexes:
                matched_donor_tx = find_matching_donor(tx, donor_index)
                if matched_donor_tx is not None:
                    matched_source = donor_basename
                    break

            if matched_donor_tx is not None:
                # Graft UTR exons onto the acceptor
                new_exons = add_utr_to_transcript(
                    acceptor_tx=tx,
                    donor_tx=matched_donor_tx,
                    max_5prime=args.max_5prime,
                    max_3prime=args.max_3prime,
                    min_exon=args.min_utr_exon,
                )
                tx.exons = new_exons
                tx.start = new_exons[0].start  if new_exons else tx.start
                tx.end   = new_exons[-1].end   if new_exons else tx.end
                # Tag transcript with its UTR source
                tx.attrs = _set_attr(tx.attrs, 'utr_source', matched_source)
                n_matched += 1
            else:
                # No donor match — keep the acceptor's existing exon structure
                # (CDS exons become the exon set if there are no exons yet)
                if not tx.exons and tx.cds:
                    tx.exons = [Exon(start=c.start, end=c.end) for c in tx.cds]

    # ------------------------------------------------------------------ #
    #  Write output                                                       #
    # ------------------------------------------------------------------ #
    g_count, tx_count = write_gff3(acceptor_genes, args.out)

    print(
        f'add_utrs: {n_acceptors} acceptor transcripts processed\n'
        f'  UTR added:   {n_matched}\n'
        f'  No CDS:      {n_skipped_no_cds}\n'
        f'  Unmatched:   {n_acceptors - n_matched - n_skipped_no_cds}\n'
        f'  Output:      {g_count} genes / {tx_count} transcripts → {args.out}',
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
