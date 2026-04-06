"""
Compute annotation statistics from a GFF3 file.

These statistics are used to compare our Python/Nextflow pipeline outputs
against reference outputs (published Ensembl GFF3, or outputs from the
Perl/eHive pipeline).

Key metrics:
  - Gene counts by biotype
  - Transcript count per gene distribution
  - UTR coverage rates (% of protein-coding genes with 5'/3' UTR)
  - Pseudogene rates
  - Readthrough rate
  - Canonical transcript CDS length distribution
  - Exon count distribution
"""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class TranscriptStats:
    biotype: str
    n_exons: int
    has_utr5: bool
    has_utr3: bool
    cds_length: int          # total CDS bp
    is_canonical: bool
    is_pseudogene: bool
    is_readthrough: bool

@dataclass
class GeneStats:
    biotype: str
    n_transcripts: int
    has_canonical: bool
    transcripts: List[TranscriptStats] = field(default_factory=list)


@dataclass
class AnnotationStats:
    """Computed statistics for one GFF3 file."""
    total_genes: int = 0
    total_transcripts: int = 0
    genes_by_biotype: Dict[str, int] = field(default_factory=Counter)
    transcripts_by_biotype: Dict[str, int] = field(default_factory=Counter)

    # UTR coverage
    protein_coding_genes: int = 0
    genes_with_utr5: int = 0       # has ≥1 transcript with 5' UTR
    genes_with_utr3: int = 0
    genes_with_both_utrs: int = 0

    # Pseudogenes
    pseudogene_count: int = 0
    processed_pseudogene_count: int = 0

    # Readthrough
    readthrough_transcript_count: int = 0

    # Transcripts per gene
    tx_per_gene: List[int] = field(default_factory=list)

    # CDS lengths of canonical transcripts
    canonical_cds_lengths: List[int] = field(default_factory=list)

    # Exon counts
    exon_counts: List[int] = field(default_factory=list)

    # Canonical coverage
    genes_with_canonical: int = 0

    def utr5_rate(self) -> float:
        if not self.protein_coding_genes:
            return 0.0
        return self.genes_with_utr5 / self.protein_coding_genes

    def utr3_rate(self) -> float:
        if not self.protein_coding_genes:
            return 0.0
        return self.genes_with_utr3 / self.protein_coding_genes

    def pseudogene_rate(self) -> float:
        if not self.total_genes:
            return 0.0
        return (self.pseudogene_count + self.processed_pseudogene_count) / self.total_genes

    def mean_tx_per_gene(self) -> float:
        return statistics.mean(self.tx_per_gene) if self.tx_per_gene else 0.0

    def median_cds_length(self) -> float:
        return statistics.median(self.canonical_cds_lengths) if self.canonical_cds_lengths else 0.0

    def canonical_coverage(self) -> float:
        if not self.total_genes:
            return 0.0
        return self.genes_with_canonical / self.total_genes

    def summary(self) -> Dict:
        return {
            "total_genes": self.total_genes,
            "total_transcripts": self.total_transcripts,
            "genes_by_biotype": dict(self.genes_by_biotype),
            "protein_coding_genes": self.protein_coding_genes,
            "utr5_rate": round(self.utr5_rate(), 4),
            "utr3_rate": round(self.utr3_rate(), 4),
            "pseudogene_rate": round(self.pseudogene_rate(), 4),
            "readthrough_transcript_count": self.readthrough_transcript_count,
            "mean_tx_per_gene": round(self.mean_tx_per_gene(), 2),
            "median_canonical_cds_bp": self.median_cds_length(),
            "canonical_coverage": round(self.canonical_coverage(), 4),
        }


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def compute_stats(gff3_path: str) -> AnnotationStats:
    """
    Parse a GFF3 file and compute AnnotationStats.

    Handles feature types: gene, mRNA/transcript, exon, CDS, UTR/five_prime_UTR,
    three_prime_UTR, pseudogenic_transcript.
    """
    stats = AnnotationStats()

    # --- Pass 1: collect features ---
    genes: Dict[str, Dict] = {}          # gene_id → {biotype, tx_ids}
    transcripts: Dict[str, Dict] = {}   # tx_id → {biotype, gene_id, exons, cds, utr5, utr3, attrs}

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            _, _, feature, start, end, _, strand, _, attr_str = cols
            start, end = int(start), int(end)
            attrs = _parse_attrs(attr_str)
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            biotype = attrs.get("biotype", "")

            if feature == "gene":
                genes[fid] = {"biotype": biotype, "tx_ids": []}

            elif feature in ("mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript",
                             "processed_pseudogene"):
                biotype_tx = attrs.get("biotype", feature)
                transcripts[fid] = {
                    "gene_id": parent, "biotype": biotype_tx,
                    "n_exons": 0, "cds_bp": 0,
                    "has_utr5": False, "has_utr3": False,
                    "is_canonical": attrs.get("canonical_transcript") == "1",
                    "is_pseudogene": biotype_tx in ("pseudogene", "processed_pseudogene",
                                                    "transcribed_pseudogene"),
                    "is_readthrough": "readthrough_transcript" in attrs,
                    "strand": strand,
                }
                if parent in genes:
                    genes[parent]["tx_ids"].append(fid)

            elif feature == "exon":
                if parent in transcripts:
                    transcripts[parent]["n_exons"] += 1

            elif feature == "CDS":
                if parent in transcripts:
                    transcripts[parent]["cds_bp"] += end - start + 1

            elif feature in ("five_prime_UTR", "five_prime_utr"):
                if parent in transcripts:
                    transcripts[parent]["has_utr5"] = True

            elif feature in ("three_prime_UTR", "three_prime_utr"):
                if parent in transcripts:
                    transcripts[parent]["has_utr3"] = True

    # --- Infer UTRs from exon > CDS if no explicit UTR features ---
    # (not computed here — requires coordinate comparison; done in comparison tests)

    # --- Pass 2: aggregate into AnnotationStats ---
    stats.total_genes = len(genes)
    stats.total_transcripts = len(transcripts)

    for g_id, gene in genes.items():
        biotype = gene["biotype"] or "unknown"
        stats.genes_by_biotype[biotype] += 1

        tx_ids = gene["tx_ids"]
        stats.tx_per_gene.append(len(tx_ids))

        is_pc = biotype == "protein_coding"
        if is_pc:
            stats.protein_coding_genes += 1

        gene_has_utr5 = gene_has_utr3 = False
        gene_has_canonical = False
        canonical_cds = 0

        for tx_id in tx_ids:
            tx = transcripts.get(tx_id)
            if not tx:
                continue
            stats.transcripts_by_biotype[tx["biotype"]] += 1
            stats.exon_counts.append(tx["n_exons"])

            if tx["is_pseudogene"]:
                if tx["biotype"] == "processed_pseudogene":
                    stats.processed_pseudogene_count += 1
                else:
                    stats.pseudogene_count += 1

            if tx["is_readthrough"]:
                stats.readthrough_transcript_count += 1

            if tx["has_utr5"] or tx["is_canonical"]:
                gene_has_utr5 = gene_has_utr5 or tx["has_utr5"]
            if tx["has_utr3"] or tx["is_canonical"]:
                gene_has_utr3 = gene_has_utr3 or tx["has_utr3"]

            if tx["is_canonical"]:
                gene_has_canonical = True
                canonical_cds = tx["cds_bp"]

        if is_pc:
            if gene_has_utr5:
                stats.genes_with_utr5 += 1
            if gene_has_utr3:
                stats.genes_with_utr3 += 1
            if gene_has_utr5 and gene_has_utr3:
                stats.genes_with_both_utrs += 1

        if gene_has_canonical:
            stats.genes_with_canonical += 1
            if canonical_cds:
                stats.canonical_cds_lengths.append(canonical_cds)

    return stats
