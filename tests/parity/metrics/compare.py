"""
GFF3 comparison utilities for biological parity testing.

Compares outputs from two annotation pipelines to find where they agree
and disagree. Used to identify whether differences between our Python/NF
implementation and the Perl/eHive reference are:
  (a) regressions — we broke something
  (b) improvements — we fixed a bug or extended an algorithm
  (c) expected differences — algorithm was intentionally changed

Key comparisons:
  - Gene overlap rate: % of test genes that have a matching reference gene
    (by CDS overlap on same strand)
  - CDS boundary exact match rate: % of matching gene pairs where all CDS
    boundaries are identical
  - UTR extension: distribution of UTR bp added vs reference
  - Biotype concordance: % of genes with same biotype in both sets
  - Pseudogene recall: % of reference pseudogenes also flagged in test
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class GffInterval:
    chrom: str
    start: int    # 1-based inclusive
    end: int
    strand: int
    feature_id: str
    biotype: str = ""


@dataclass
class ComparisonResult:
    """Results of comparing a test GFF3 against a reference GFF3."""

    # Gene-level
    total_ref_genes: int = 0
    total_test_genes: int = 0
    matched_genes: int = 0          # ref genes with a test gene overlap (CDS)
    unmatched_ref_genes: int = 0    # ref genes with no test equivalent
    novel_test_genes: int = 0       # test genes with no ref equivalent

    # CDS boundary accuracy (among matched genes)
    exact_cds_match: int = 0        # test CDS boundaries == ref CDS boundaries
    partial_cds_match: int = 0      # some but not all CDS boundaries match
    cds_boundary_accuracy: float = 0.0  # exact_cds / matched

    # Biotype concordance
    biotype_concordant: int = 0
    biotype_discordant: int = 0

    # UTR changes (test vs reference)
    utr5_added: List[int] = field(default_factory=list)   # bp of 5' UTR added (positive = more in test)
    utr3_added: List[int] = field(default_factory=list)

    # Pseudogene metrics
    ref_pseudogenes: int = 0
    test_pseudogenes: int = 0
    pseudogene_recall: float = 0.0   # test pseudo / ref pseudo (among matched)
    pseudogene_precision: float = 0.0  # ref pseudo / test pseudo

    def summary(self) -> Dict:
        gene_recall = self.matched_genes / self.total_ref_genes if self.total_ref_genes else 0.0
        gene_precision = self.matched_genes / self.total_test_genes if self.total_test_genes else 0.0
        return {
            "ref_genes": self.total_ref_genes,
            "test_genes": self.total_test_genes,
            "matched_genes": self.matched_genes,
            "gene_recall": round(gene_recall, 4),
            "gene_precision": round(gene_precision, 4),
            "unmatched_ref_genes": self.unmatched_ref_genes,
            "novel_test_genes": self.novel_test_genes,
            "cds_boundary_accuracy": round(self.cds_boundary_accuracy, 4),
            "biotype_concordance": round(
                self.biotype_concordant / max(self.matched_genes, 1), 4),
            "mean_utr5_delta_bp": (
                round(sum(self.utr5_added) / len(self.utr5_added), 1)
                if self.utr5_added else None
            ),
            "mean_utr3_delta_bp": (
                round(sum(self.utr3_added) / len(self.utr3_added), 1)
                if self.utr3_added else None
            ),
            "pseudogene_recall": round(self.pseudogene_recall, 4),
            "pseudogene_precision": round(self.pseudogene_precision, 4),
        }


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


def _load_genes(gff3_path: str) -> Dict[str, Dict]:
    """
    Load genes from GFF3 into a dict: gene_id → {chrom, start, end, strand,
    biotype, cds_intervals: [(start, end)], utr5_bp, utr3_bp, is_pseudogene}.

    CDS intervals are the union of all CDS features across all transcripts
    (used for gene-level overlap detection).
    """
    genes: Dict[str, Dict] = {}
    transcripts: Dict[str, str] = {}  # tx_id → gene_id

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr_str = cols
            start, end = int(start), int(end)
            attrs = _parse_attrs(attr_str)
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            biotype = attrs.get("biotype", "")

            if feature == "gene":
                genes[fid] = {
                    "chrom": chrom, "start": start, "end": end,
                    "strand": 1 if strand == "+" else -1,
                    "biotype": biotype, "cds_intervals": [],
                    "utr5_bp": 0, "utr3_bp": 0,
                    "is_pseudogene": biotype in ("pseudogene", "processed_pseudogene",
                                                  "transcribed_pseudogene"),
                    "tx_ids": [],
                }

            elif feature in ("mRNA", "transcript", "lnc_RNA", "pseudogenic_transcript"):
                transcripts[fid] = parent
                if parent in genes:
                    genes[parent]["tx_ids"].append(fid)

            elif feature == "CDS":
                gene_id = transcripts.get(parent)
                if gene_id and gene_id in genes:
                    genes[gene_id]["cds_intervals"].append((start, end))

            elif feature in ("five_prime_UTR", "five_prime_utr"):
                gene_id = transcripts.get(parent)
                if gene_id and gene_id in genes:
                    genes[gene_id]["utr5_bp"] += end - start + 1

            elif feature in ("three_prime_UTR", "three_prime_utr"):
                gene_id = transcripts.get(parent)
                if gene_id and gene_id in genes:
                    genes[gene_id]["utr3_bp"] += end - start + 1

    return genes


def _genes_overlap(g1: Dict, g2: Dict) -> bool:
    """
    Two genes overlap if they share the same chrom+strand and their
    CDS intervals have any genomic overlap.
    """
    if g1["chrom"] != g2["chrom"] or g1["strand"] != g2["strand"]:
        return False
    for (s1, e1) in g1["cds_intervals"]:
        for (s2, e2) in g2["cds_intervals"]:
            if s1 <= e2 and s2 <= e1:
                return True
    return False


def _cds_exact_match(g1: Dict, g2: Dict) -> bool:
    """True if the sorted CDS interval sets are identical."""
    set1 = frozenset(g1["cds_intervals"])
    set2 = frozenset(g2["cds_intervals"])
    return set1 == set2


def compare_gff3(ref_path: str, test_path: str) -> ComparisonResult:
    """
    Compare a test GFF3 against a reference GFF3.

    Matching is done by CDS overlap (not by gene name) to handle cases
    where gene IDs differ between the two pipelines.

    Returns a ComparisonResult with gene-level and transcript-level metrics.
    """
    result = ComparisonResult()

    ref_genes = _load_genes(ref_path)
    test_genes = _load_genes(test_path)

    result.total_ref_genes = len(ref_genes)
    result.total_test_genes = len(test_genes)
    result.ref_pseudogenes = sum(1 for g in ref_genes.values() if g["is_pseudogene"])
    result.test_pseudogenes = sum(1 for g in test_genes.values() if g["is_pseudogene"])

    # Build a spatial index: chrom+strand → list of test genes (simple linear scan
    # is fine for parity tests on small genomes; use interval tree for large ones)
    test_by_chrom: Dict[Tuple, List] = defaultdict(list)
    for gid, g in test_genes.items():
        test_by_chrom[(g["chrom"], g["strand"])].append((gid, g))

    matched_test_ids: Set[str] = set()
    matched_ref_pseudogenes = 0
    matched_test_pseudogenes = 0

    for ref_id, ref_gene in ref_genes.items():
        candidates = test_by_chrom.get((ref_gene["chrom"], ref_gene["strand"]), [])
        match = None
        for test_id, test_gene in candidates:
            if _genes_overlap(ref_gene, test_gene):
                match = (test_id, test_gene)
                break

        if match is None:
            result.unmatched_ref_genes += 1
            continue

        test_id, test_gene = match
        result.matched_genes += 1
        matched_test_ids.add(test_id)

        # CDS boundary accuracy
        if _cds_exact_match(ref_gene, test_gene):
            result.exact_cds_match += 1
        else:
            result.partial_cds_match += 1

        # Biotype concordance
        if ref_gene["biotype"] == test_gene["biotype"]:
            result.biotype_concordant += 1
        else:
            result.biotype_discordant += 1

        # UTR deltas (positive = test has more UTR)
        result.utr5_added.append(test_gene["utr5_bp"] - ref_gene["utr5_bp"])
        result.utr3_added.append(test_gene["utr3_bp"] - ref_gene["utr3_bp"])

        # Pseudogene concordance
        if ref_gene["is_pseudogene"]:
            matched_ref_pseudogenes += 1
            if test_gene["is_pseudogene"]:
                matched_test_pseudogenes += 1

    result.novel_test_genes = len(test_genes) - len(matched_test_ids)
    result.cds_boundary_accuracy = (
        result.exact_cds_match / result.matched_genes if result.matched_genes else 0.0
    )
    result.pseudogene_recall = (
        matched_test_pseudogenes / matched_ref_pseudogenes
        if matched_ref_pseudogenes else 0.0
    )
    result.pseudogene_precision = (
        matched_ref_pseudogenes / result.test_pseudogenes
        if result.test_pseudogenes else 0.0
    )

    return result
