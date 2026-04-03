#!/usr/bin/env python3
"""
collapse_long_reads.py — Collapse overlapping long-read BAM alignments
into non-redundant GFF3 transcript models.

Mirrors the logic of HiveTranscriptCoalescer but operates entirely on
flat files. Outputs GFF3 with biotype=isoseq (ONT/cDNA) or biotype=cdna.

Usage:
    collapse_long_reads.py --bam in.bam --out out.gff3 --stats out.stats.tsv
                           [--min-overlap 0.8] [--max-intron 200000]
                           [--sample-id SAMPLE]
"""

import argparse
import sys
import collections
from pathlib import Path

try:
    import pysam
except ImportError:
    sys.exit("pysam is required: pip install pysam")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class Exon:
    __slots__ = ('start', 'end')

    def __init__(self, start: int, end: int):
        self.start = start
        self.end   = end

    def length(self) -> int:
        return self.end - self.start

    def overlap(self, other: 'Exon') -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))


class Transcript:
    def __init__(self, tid: str, chrom: str, strand: str, exons: list[Exon]):
        self.tid    = tid
        self.chrom  = chrom
        self.strand = strand
        self.exons  = sorted(exons, key=lambda e: e.start)

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end

    def exon_length(self) -> int:
        return sum(e.length() for e in self.exons)

    def reciprocal_overlap(self, other: 'Transcript') -> float:
        """Reciprocal exon-level overlap fraction (0–1)."""
        if self.chrom != other.chrom or self.strand != other.strand:
            return 0.0
        shared = 0
        for e1 in self.exons:
            for e2 in other.exons:
                shared += e1.overlap(e2)
        total = min(self.exon_length(), other.exon_length())
        return shared / total if total > 0 else 0.0


# ---------------------------------------------------------------------------
# BAM parsing
# ---------------------------------------------------------------------------

def parse_bam(bam_path: str, max_intron: int) -> list[Transcript]:
    """Extract splice-aware transcript models from a BAM file."""
    transcripts = []
    seen_ids    = set()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for i, read in enumerate(bam.fetch()):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            chrom  = bam.get_reference_name(read.reference_id)
            strand = '-' if read.is_reverse else '+'
            exons  = _cigar_to_exons(read.reference_start, read.cigartuples or [], max_intron)

            if not exons:
                continue

            tid = f"{chrom}_{strand}_{read.reference_start}_{read.reference_end}_{i}"
            if tid in seen_ids:
                continue
            seen_ids.add(tid)

            transcripts.append(Transcript(tid, chrom, strand, exons))

    return transcripts


def _cigar_to_exons(ref_start: int, cigar: list, max_intron: int) -> list[Exon]:
    """Convert a CIGAR string to a list of exon intervals (reference coords)."""
    # CIGAR op codes: 0=M, 1=I, 2=D, 3=N(intron), 4=S, 5=H, 6=P, 7==, 8=X
    CONSUMES_REF = {0, 2, 3, 7, 8}
    INTRON_OP    = 3

    exons   = []
    pos     = ref_start
    ex_start = ref_start
    in_exon  = True

    for op, length in cigar:
        if op == INTRON_OP:
            if in_exon and pos > ex_start:
                exons.append(Exon(ex_start, pos))
            if length > max_intron:
                return []  # filter oversized introns
            pos     += length
            ex_start = pos
            in_exon  = True
        elif op in CONSUMES_REF:
            pos += length

    if in_exon and pos > ex_start:
        exons.append(Exon(ex_start, pos))

    return exons if exons else []


# ---------------------------------------------------------------------------
# Clustering and collapse
# ---------------------------------------------------------------------------

def cluster_by_locus(transcripts: list[Transcript]) -> list[list[Transcript]]:
    """Group transcripts into overlapping loci (per chrom+strand)."""
    by_strand: dict[tuple, list[Transcript]] = collections.defaultdict(list)
    for t in transcripts:
        by_strand[(t.chrom, t.strand)].append(t)

    clusters = []
    for key, txs in by_strand.items():
        txs_sorted = sorted(txs, key=lambda t: (t.start, t.end))
        current    = [txs_sorted[0]]
        max_end    = txs_sorted[0].end

        for tx in txs_sorted[1:]:
            if tx.start < max_end:  # overlaps current cluster
                current.append(tx)
                max_end = max(max_end, tx.end)
            else:
                clusters.append(current)
                current = [tx]
                max_end = tx.end

        clusters.append(current)

    return clusters


def collapse_cluster(cluster: list[Transcript], min_overlap: float) -> list[Transcript]:
    """Greedily collapse transcripts within a cluster by reciprocal exon overlap."""
    if not cluster:
        return []

    # Sort longest first — use the longest as the representative
    sorted_txs = sorted(cluster, key=lambda t: t.exon_length(), reverse=True)
    representatives: list[Transcript] = []

    for tx in sorted_txs:
        merged = False
        for rep in representatives:
            if tx.reciprocal_overlap(rep) >= min_overlap:
                merged = True
                break
        if not merged:
            representatives.append(tx)

    return representatives


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------

def write_gff3(transcripts: list[Transcript], out_path: str, sample_id: str) -> None:
    gene_counter = 1
    # Group transcripts into genes by locus (overlapping same-strand models)
    clusters = cluster_by_locus(transcripts)

    with open(out_path, 'w') as fh:
        fh.write("##gff-version 3\n")
        for cluster in clusters:
            if not cluster:
                continue
            chrom  = cluster[0].chrom
            strand = cluster[0].strand
            g_start = min(t.start for t in cluster) + 1  # GFF3 is 1-based
            g_end   = max(t.end   for t in cluster)
            gene_id = f"{sample_id}_gene_{gene_counter:06d}"
            gene_counter += 1

            fh.write(
                f"{chrom}\tlong_read\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t"
                f"ID={gene_id};biotype=isoseq\n"
            )

            for tx_i, tx in enumerate(cluster, start=1):
                tx_start = tx.start + 1
                tx_end   = tx.end
                tx_id    = f"{gene_id}.t{tx_i}"
                fh.write(
                    f"{chrom}\tlong_read\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t"
                    f"ID={tx_id};Parent={gene_id};biotype=isoseq\n"
                )
                for ex in tx.exons:
                    fh.write(
                        f"{chrom}\tlong_read\texon\t{ex.start + 1}\t{ex.end}\t.\t{strand}\t.\t"
                        f"Parent={tx_id}\n"
                    )


def write_stats(transcripts: list[Transcript], out_path: str, sample_id: str) -> None:
    clusters = cluster_by_locus(transcripts)
    n_genes  = len(clusters)
    n_tx     = len(transcripts)

    with open(out_path, 'w') as fh:
        fh.write("sample\tgenes\ttranscripts\n")
        fh.write(f"{sample_id}\t{n_genes}\t{n_tx}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam',         required=True,  help='Input BAM file')
    parser.add_argument('--out',         required=True,  help='Output GFF3 file')
    parser.add_argument('--stats',       required=True,  help='Output stats TSV')
    parser.add_argument('--min-overlap', type=float, default=0.8,
                        help='Min reciprocal exon overlap to merge transcripts (default: 0.8)')
    parser.add_argument('--max-intron',  type=int,   default=200000,
                        help='Max intron size in bp; larger alignments are filtered (default: 200000)')
    parser.add_argument('--sample-id',   default='sample',
                        help='Sample identifier used in feature IDs (default: sample)')
    args = parser.parse_args()

    print(f"[collapse_long_reads] Parsing BAM: {args.bam}", flush=True)
    transcripts = parse_bam(args.bam, args.max_intron)
    print(f"[collapse_long_reads] Parsed {len(transcripts)} alignments", flush=True)

    clusters  = cluster_by_locus(transcripts)
    collapsed = []
    for cluster in clusters:
        collapsed.extend(collapse_cluster(cluster, args.min_overlap))

    print(f"[collapse_long_reads] Collapsed to {len(collapsed)} transcript models "
          f"from {len(clusters)} loci", flush=True)

    write_gff3(collapsed, args.out, args.sample_id)
    write_stats(collapsed, args.stats, args.sample_id)
    print(f"[collapse_long_reads] Written to {args.out}", flush=True)


if __name__ == '__main__':
    main()
