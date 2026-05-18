#!/usr/bin/env python3
"""
Build a translon_db manifest TSV from raw tool outputs.

Handles all 5 tools (PRICE, RiboTIE, ORFQuant, iRibo, RibORF2) with
CDS / non-CDS split from each tool's native output structure.

The manifest is passed directly to translon_db_standardise.py via --manifest,
so original file paths are preserved and source_feature_class is inferred from
the original filenames (annotated/known -> cds, novel/unannotated -> non_cds).

Usage:
    python build_translon_manifest.py \\
        --results-dir /hps/.../full_pilot_results \\
        --riborf2-converted /hps/.../riborf2_converted \\
        --out manifest.tsv

    # Run DB build directly against original files:
    translon_db_standardise.py \\
        --input-root /any/valid/dir \\
        --manifest manifest.tsv \\
        --out-dir translon_db/ \\
        --fasta genome.fa \\
        --gtf gencode.gtf.gz
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Sample name normalisation
# ---------------------------------------------------------------------------

PRICE_SAMPLE_MAP: dict[str, str] = {
    "Fibo":     "Ribo_Fib_pooled",
    "Fibo_0":   "SRR15513179",
    "Fibo_1":   "SRR15513180",
    "Fibo_2":   "SRR15513181",
    "Fibo_3":   "SRR15513182",
    "Fibo_4":   "Fib_24_45m",
    "Fibo_5":   "Fib_24_bsl",
    "Fibo_6":   "Fib_27_45m",
    "Fibo_7":   "Fib_27_bsl",
    "Fibo_8":   "Fib_41_45m",
    "Fibo_9":   "Fib_41_bsl",
    "Endo":     "Ribo_EC_pooled",
    "Endo_0":   "SRR15513197",
    "Endo_1":   "SRR15513198_GENELAB1026",
    "Endo_2":   "SRR15513199",
    "Endo_3":   "SRR15513200",
    "Endo_4":   "SRR15513201",
    "Endo_5":   "SRR15513202_GENELAB1143",
    "Pancreas":   "Ribo_Pancreas_pooled",
    "Pancreas_0": "SRR11005875_to_79",
    "Pancreas_1": "SRR11005880_to_84",
    "Pancreas_2": "SRR11005885_to_89",
    "Pancreas_3": "SRR11005890_to_94",
    "Pancreas_4": "SRR11005895_to_99",
    "Pancreas_5": "SRR11005900_to_04",
}

# iRibo uses abbreviated pooled names
IRIBO_POOLED_MAP: dict[str, str] = {
    "Ribo_EC_p":       "Ribo_EC_pooled",
    "Ribo_Fib_p":      "Ribo_Fib_pooled",
    "Ribo_pancreas_p": "Ribo_Pancreas_pooled",
}

# Capitalisation / punctuation discrepancies across tools
MISC_FIXES: dict[str, str] = {
    "Ribo_pancreas_pooled": "Ribo_Pancreas_pooled",
    "Ribo_ECs_pooled":      "Ribo_EC_pooled",
    "Ribo_fib_pooled":      "Ribo_Fib_pooled",
}

EXPECTED_SAMPLES: frozenset[str] = frozenset({
    # Fibroblast individual (GENELAB dataset — named by condition)
    "Fib_24_45m", "Fib_24_bsl",
    "Fib_27_45m", "Fib_27_bsl",
    "Fib_41_45m", "Fib_41_bsl",
    # Fibroblast individual (second dataset — named by SRR accession, PRICE Fibo_0–3)
    "SRR15513179", "SRR15513180", "SRR15513181", "SRR15513182",
    # Fibroblast pooled
    "Ribo_Fib_pooled",
    # Endothelial individual
    "SRR15513197", "SRR15513198_GENELAB1026",
    "SRR15513199", "SRR15513200",
    "SRR15513201", "SRR15513202_GENELAB1143",
    # Endothelial pooled
    "Ribo_EC_pooled",
    # Pancreas individual
    "SRR11005875_to_79", "SRR11005880_to_84", "SRR11005885_to_89",
    "SRR11005890_to_94", "SRR11005895_to_99", "SRR11005900_to_04",
    # Pancreas pooled
    "Ribo_Pancreas_pooled",
})


def normalise(raw: str) -> str:
    """Normalise a raw sample token to a canonical sample_id."""
    s = raw.strip()
    s = s.replace("GENELAB-000", "GENELAB").replace("GENELAB_000", "GENELAB")
    s = re.sub(r"_1$", "", s)          # strip trailing _1 from SRR15513179_1 etc.
    s = PRICE_SAMPLE_MAP.get(s, s)
    s = IRIBO_POOLED_MAP.get(s, s)
    s = MISC_FIXES.get(s, s)
    return s


# ---------------------------------------------------------------------------
# Per-tool file collection
# ---------------------------------------------------------------------------

def collect_price(results_dir: Path) -> list[dict]:
    """
    PRICE_results/price/{sample}.known.bed  -> CDS
    PRICE_results/price/{sample}.novel.bed  -> non-CDS
    Skip pancreas_mymapping_* (FASTQ re-run duplicates).
    """
    rows = []
    for f in sorted((results_dir / "PRICE_results" / "price").glob("*.bed")):
        m = re.match(r"^(.+)\.(known|novel)\.bed$", f.name)
        if not m:
            continue
        raw_sample = m.group(1)
        if "mymapping" in raw_sample:
            continue
        rows.append({"path": str(f), "tool": "PRICE", "sample_id": normalise(raw_sample)})
    return rows


def collect_ribotie(results_dir: Path) -> list[dict]:
    """
    RiboTIE_results/deliverables/annotated/db_*.annotated.out.gtf  -> CDS
    RiboTIE_results/deliverables/novel/db_*.novel.out.gtf           -> non-CDS

    GTF files only (Pancreas_N entries in deliverables are FASTQ-run CSV — skip).
    The unfiltered/ subdirectory is also skipped.
    """
    rows = []
    for subdir in ("annotated", "novel"):
        d = results_dir / "RiboTIE_results" / "deliverables" / subdir
        for f in sorted(d.glob("*.gtf")):
            if "unfiltered" in str(f):
                continue
            stem = f.name
            stem = re.sub(r"^db_", "", stem)
            stem = re.sub(r"\.(annotated|novel)\.out\.gtf$", "", stem)
            stem = re.sub(r"\.Aligned\.toTranscriptome\.out$", "", stem)
            stem = re.sub(r"_Transcriptome$", "", stem)   # pooled: Ribo_*_pooled_Transcriptome
            rows.append({"path": str(f), "tool": "RiboTIE", "sample_id": normalise(stem)})
    return rows


def collect_orfquant(results_dir: Path) -> list[dict]:
    """
    ORFQuant_results/bed_files/annotated_orf_exon_genomic_*.bed -> CDS
    ORFQuant_results/bed_files/novel_orf_exon_genomic_*.bed     -> non-CDS

    Some samples have two files for the same sample_id+class (one with
    .Aligned.sortedByCoord.out in the name, one without). Keep the shorter
    (simpler) filename.
    """
    bed_dir = results_dir / "ORFQuant_results" / "bed_files"
    # (sample_id, cls) -> (name_length, path)
    best: dict[tuple[str, str], tuple[int, Path]] = {}
    for f in sorted(bed_dir.glob("*.bed")):
        m = re.match(
            r"^(annotated|novel)_orf_exon_genomic_"
            r"(.+?)(?:_S\d+_R1_001(?:_trimmed)?)?(?:\.Aligned\.sortedByCoord\.out)?\.bed$",
            f.name,
        )
        if not m:
            continue
        cls, raw_sample = m.group(1), m.group(2)
        sample_id = normalise(raw_sample)
        key = (sample_id, cls)
        if key not in best or len(f.name) < best[key][0]:
            best[key] = (len(f.name), f)
    return [
        {"path": str(path), "tool": "ORFQuant", "sample_id": sid}
        for (sid, _), (_, path) in sorted(best.items())
    ]


def collect_iribo(results_dir: Path) -> list[dict]:
    """
    iRibo_results/annotated_orfs_*.bed   -> CDS
    iRibo_results/unannotated_orfs_*.bed -> non-CDS

    Fibroblast files carry _S##_R1_001 sequencing run suffix — strip it.
    Pooled samples use abbreviated names (Ribo_EC_p etc.) — normalise.
    Skip *_fastq variants.
    """
    rows = []
    for f in sorted((results_dir / "iRibo_results").glob("*.bed")):
        m = re.match(r"^(annotated|unannotated)_orfs_(.+?)(?:_S\d+_R1_001)?\.bed$", f.name)
        if not m:
            continue
        raw_sample = m.group(2)
        if raw_sample.endswith("_fastq"):
            continue
        rows.append({"path": str(f), "tool": "iRibo", "sample_id": normalise(raw_sample)})
    return rows


def collect_riborf2(converted_dir: Path) -> list[dict]:
    """
    RibORF2 outputs only repre.valid.ORF.genepred.txt (no CDS/non-CDS split).
    These must be pre-converted to BED12 via genePredToBed.
    source_feature_class will be 'unknown'; the DB reference_cds JOIN classifies them.
    """
    if not converted_dir.exists():
        print(
            f"[WARN] RibORF2 converted dir not found: {converted_dir}\n"
            f"       Run genePredToBed on RibORF_results/*/repre.valid.ORF.genepred.txt first.",
            file=sys.stderr,
        )
        return []
    rows = []
    for f in sorted(converted_dir.glob("*.bed12")):
        # Strip trailing _S##_R1_001_trimmed artefact if present
        stem = re.sub(r"_S\d+_R1_001(?:_trimmed)?$", "", f.stem)
        rows.append({"path": str(f), "tool": "RibORF2", "sample_id": normalise(stem)})
    return rows


# ---------------------------------------------------------------------------
# Summary / validation
# ---------------------------------------------------------------------------

def _classify(path_str: str) -> str:
    lower = path_str.lower()
    # Check unannotated BEFORE annotated: "unannotated" contains the substring "annotated"
    if any(x in lower for x in ("novel", "unannotated")):
        return "non_cds"
    if any(x in lower for x in ("known", "annotated")):
        return "cds"
    return "unknown"


def print_summary(rows: list[dict], expected: frozenset[str]) -> None:
    TOOLS = ["PRICE", "RiboTIE", "ORFQuant", "iRibo", "RibORF2"]
    W = 60

    print("\n" + "=" * W)
    print("MANIFEST SUMMARY")
    print("=" * W)

    # Index rows by tool
    by_tool: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_tool[r["tool"]].append(r)

    all_ok = True

    for tool in TOOLS:
        tool_rows = by_tool.get(tool, [])
        if not tool_rows:
            print(f"\n{tool}  [NO FILES FOUND]")
            all_ok = False
            continue

        by_cls: dict[str, set[str]] = defaultdict(set)
        for r in tool_rows:
            by_cls[_classify(r["path"])].add(r["sample_id"])

        cds     = by_cls["cds"]
        non_cds = by_cls["non_cds"]
        unk     = by_cls["unknown"]
        all_s   = cds | non_cds | unk

        print(f"\n{tool}  ({len(tool_rows)} files, {len(all_s)} samples)")
        print(f"  CDS files     : {len(cds):2d} samples")
        print(f"  non-CDS files : {len(non_cds):2d} samples")
        if unk:
            print(f"  unknown class : {len(unk):2d} samples  (no keyword in filename)")

        missing = expected - all_s
        extra   = all_s - expected
        if missing:
            print(f"  MISSING from expected set ({len(missing)}): {sorted(missing)}")
            all_ok = False
        if extra:
            print(f"  Extra / unexpected ({len(extra)}): {sorted(extra)}")

        if cds and non_cds:
            cds_only     = sorted(cds - non_cds)
            non_cds_only = sorted(non_cds - cds)
            if cds_only:
                print(f"  CDS only (no non-CDS pair)    : {cds_only}")
            if non_cds_only:
                print(f"  non-CDS only (no CDS pair)    : {non_cds_only}")

        if not missing and not extra:
            print(f"  Coverage: complete ({len(all_s)}/{len(expected)} expected samples)")

    # Cross-tool matrix
    print(f"\n{'─' * W}")
    print("CROSS-TOOL SAMPLE COVERAGE")
    print(f"{'Sample':<30}", end="")
    for t in TOOLS:
        print(f"  {t[:8]:<8}", end="")
    print()
    print("─" * W)

    samples_per_tool = {t: {r["sample_id"] for r in by_tool.get(t, [])} for t in TOOLS}
    for sample in sorted(expected):
        print(f"  {sample:<28}", end="")
        for t in TOOLS:
            present = sample in samples_per_tool[t]
            print(f"  {'✓' if present else '✗':<8}", end="")
        print()

    # Totals row
    print("─" * W)
    print(f"  {'TOTAL':<28}", end="")
    for t in TOOLS:
        n = len(samples_per_tool[t] & expected)
        print(f"  {n}/{len(expected):<6}", end="")
    print()

    print(f"\nTotal manifest entries : {len(rows)}")
    print(f"Overall status         : {'OK' if all_ok else 'INCOMPLETE — see warnings above'}")
    print("=" * W)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, required=True,
                   help="Path to full_pilot_results directory")
    p.add_argument("--riborf2-converted", type=Path, required=True,
                   help="Directory of RibORF2 BED12 files converted from GenePred")
    p.add_argument("--out", type=Path, required=True,
                   help="Output manifest TSV (path / tool / sample_id)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    collectors = [
        ("PRICE",    collect_price,    (args.results_dir,)),
        ("RiboTIE",  collect_ribotie,  (args.results_dir,)),
        ("ORFQuant", collect_orfquant, (args.results_dir,)),
        ("iRibo",    collect_iribo,    (args.results_dir,)),
        ("RibORF2",  collect_riborf2,  (args.riborf2_converted,)),
    ]

    rows: list[dict] = []
    for label, fn, fn_args in collectors:
        tool_rows = fn(*fn_args)
        rows.extend(tool_rows)
        print(f"  {label:<10}: {len(tool_rows):3d} files")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["path", "tool", "sample_id"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} entries -> {args.out}")
    print_summary(rows, EXPECTED_SAMPLES)


if __name__ == "__main__":
    main()
