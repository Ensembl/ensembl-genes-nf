#!/usr/bin/env python3
"""
Run a biological parity comparison between a reference GFF3 and our new pipeline output.

Wraps tests/parity/metrics/compare.py and gff3_stats.py into a single CLI tool
that produces a human-readable report and a JSON results file.

Typical usage:

  # After downloading the Ensembl reference GFF3:
  python run_parity_comparison.py \\
      --ref reference.gff3 \\
      --test /hps/scratch/.../finalise_geneset/finalise_geneset/final.canonical.gff3 \\
      --output parity_results_GCA_003957565.2.json \\
      --label "Nextflow v1 vs Ensembl release 113"

  # Multiple test GFF3s (e.g. to compare intermediate stages):
  python run_parity_comparison.py \\
      --ref reference.gff3 \\
      --test utr_added.gff3 finalise.gff3 \\
      --labels "post-UTR" "post-finalise" \\
      --output parity_results.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow running from the repo root or scripts/ dir
_repo_root = Path(__file__).parent.parent
sys.path.insert(0, str(_repo_root / "tests" / "parity"))
sys.path.insert(0, str(_repo_root / "tests"))

from metrics.compare import compare_gff3, ComparisonResult
from metrics.gff3_stats import compute_stats, AnnotationStats


def _format_pct(v: float) -> str:
    return f"{v * 100:.1f}%"


def _report_stats(stats: AnnotationStats, label: str) -> None:
    print(f"\n  {label}:")
    print(f"    Genes:                {stats.total_genes:,}")
    print(f"    Transcripts:          {stats.total_transcripts:,}")
    print(f"    Protein-coding genes: {stats.protein_coding_genes:,}")
    print(f"    UTR5 coverage:        {_format_pct(stats.utr5_rate())}")
    print(f"    UTR3 coverage:        {_format_pct(stats.utr3_rate())}")
    print(f"    Pseudogene rate:      {_format_pct(stats.pseudogene_rate())}")
    print(f"    Mean tx/gene:         {stats.mean_tx_per_gene():.2f}")
    print(f"    Median canonical CDS: {stats.median_cds_length():,.0f} bp")
    print(f"    Canonical coverage:   {_format_pct(stats.canonical_coverage())}")


def _report_comparison(result: ComparisonResult, label: str) -> None:
    s = result.summary()
    print(f"\n  Comparison ({label}):")
    print(f"    Gene recall:          {_format_pct(s['gene_recall'])}  "
          f"({result.matched_genes}/{result.total_ref_genes} ref genes matched)")
    print(f"    Gene precision:       {_format_pct(s['gene_precision'])}  "
          f"({result.matched_genes}/{result.total_test_genes} test genes matched)")
    print(f"    Unmatched ref genes:  {result.unmatched_ref_genes:,}  (missed)")
    print(f"    Novel test genes:     {result.novel_test_genes:,}  (new calls)")
    print(f"    CDS boundary accuracy:{_format_pct(s['cds_boundary_accuracy'])}  "
          f"(exact CDS match among matched genes)")
    print(f"    Biotype concordance:  {_format_pct(s['biotype_concordance'])}")
    if s["mean_utr5_delta_bp"] is not None:
        sign5 = "+" if s["mean_utr5_delta_bp"] >= 0 else ""
        sign3 = "+" if s["mean_utr3_delta_bp"] >= 0 else ""
        print(f"    Mean 5' UTR delta:    {sign5}{s['mean_utr5_delta_bp']:.0f} bp  "
              f"(positive = test has more)")
        print(f"    Mean 3' UTR delta:    {sign3}{s['mean_utr3_delta_bp']:.0f} bp")
    print(f"    Pseudogene recall:    {_format_pct(s['pseudogene_recall'])}")
    print(f"    Pseudogene precision: {_format_pct(s['pseudogene_precision'])}")


def _classify_differences(result: ComparisonResult) -> list[str]:
    """
    Return a list of natural-language statements about what the differences mean.
    This helps distinguish regressions from improvements from intentional changes.
    """
    notes = []
    s = result.summary()

    if s["gene_recall"] >= 0.99:
        notes.append("✓ Gene recall near-perfect — no significant missed genes")
    elif s["gene_recall"] >= 0.95:
        notes.append(f"△ Gene recall {_format_pct(s['gene_recall'])} — small number of missed genes, investigate")
    else:
        notes.append(f"✗ Gene recall {_format_pct(s['gene_recall'])} — significant gene loss, likely regression")

    if s["gene_precision"] < 0.98:
        notes.append(
            f"△ Gene precision {_format_pct(s['gene_precision'])} — "
            f"{result.novel_test_genes} novel calls. Could be improvements (novel genes) "
            "or spurious calls — check loci manually."
        )

    if s["cds_boundary_accuracy"] >= 0.99:
        notes.append("✓ CDS boundaries near-identical — algorithm parity confirmed")
    elif s["cds_boundary_accuracy"] >= 0.95:
        notes.append(
            f"△ CDS boundary accuracy {_format_pct(s['cds_boundary_accuracy'])} — "
            "small differences. Check if these correspond to known Perl bugs or new regressions."
        )
    else:
        notes.append(
            f"✗ CDS boundary accuracy {_format_pct(s['cds_boundary_accuracy'])} — "
            "major CDS boundary differences. Likely algorithm regression."
        )

    utr5 = s["mean_utr5_delta_bp"]
    utr3 = s["mean_utr3_delta_bp"]
    if utr5 is not None:
        if utr5 > 50:
            notes.append(
                f"↑ Mean +{utr5:.0f} bp 5' UTR, +{utr3:.0f} bp 3' UTR — "
                "new pipeline adds more UTR. Expected if utr_addition is improved."
            )
        elif utr5 < -50:
            notes.append(
                f"↓ Mean {utr5:.0f} bp 5' UTR — "
                "new pipeline has LESS UTR than reference. Possible regression."
            )
        else:
            notes.append(f"✓ UTR lengths similar (delta <50 bp on average)")

    if s["pseudogene_recall"] < 0.9 and result.ref_pseudogenes > 0:
        notes.append(
            f"✗ Pseudogene recall {_format_pct(s['pseudogene_recall'])} — "
            "new pipeline misses pseudogenes. Regression in pseudogene filter."
        )
    elif result.ref_pseudogenes > 0:
        notes.append(f"✓ Pseudogene recall {_format_pct(s['pseudogene_recall'])}")

    return notes


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Biological parity comparison between reference and test GFF3.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--ref", required=True, help="Reference GFF3 (e.g. published Ensembl GFF3)")
    parser.add_argument("--test", nargs="+", required=True,
                        help="One or more test GFF3 files (new pipeline outputs)")
    parser.add_argument("--labels", nargs="*",
                        help="Labels for each test GFF3 (defaults to filename)")
    parser.add_argument("--output", default="parity_results.json",
                        help="JSON output path for detailed results")
    parser.add_argument("--label", default="",
                        help="Overall experiment label for the report")

    args = parser.parse_args()

    labels = args.labels or [Path(t).name for t in args.test]
    if len(labels) < len(args.test):
        labels += [Path(t).name for t in args.test[len(labels):]]

    print("=" * 65)
    print(f"  Biological Parity Comparison")
    if args.label:
        print(f"  {args.label}")
    print("=" * 65)

    # Reference stats
    print(f"\nReference GFF3: {args.ref}")
    try:
        ref_stats = compute_stats(args.ref)
        _report_stats(ref_stats, "Reference annotation stats")
        ref_stats_dict = ref_stats.summary()
    except Exception as e:
        print(f"  WARNING: Could not compute ref stats: {e}")
        ref_stats_dict = {}

    all_results = {
        "experiment_label": args.label,
        "reference": args.ref,
        "reference_stats": ref_stats_dict,
        "comparisons": [],
    }

    for test_path, label in zip(args.test, labels):
        print(f"\nTest GFF3: {test_path}  [{label}]")
        try:
            test_stats = compute_stats(test_path)
            _report_stats(test_stats, f"Test annotation stats [{label}]")
            test_stats_dict = test_stats.summary()
        except Exception as e:
            print(f"  WARNING: Could not compute test stats: {e}")
            test_stats_dict = {}

        try:
            result = compare_gff3(args.ref, test_path)
            _report_comparison(result, label)

            notes = _classify_differences(result)
            print(f"\n  Interpretation ({label}):")
            for note in notes:
                print(f"    {note}")

            all_results["comparisons"].append({
                "label": label,
                "test_gff3": test_path,
                "test_stats": test_stats_dict,
                "comparison": result.summary(),
                "interpretation": notes,
            })

        except Exception as e:
            print(f"  ERROR during comparison: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'=' * 65}")
    print(f"  Results written to: {args.output}")
    print("=" * 65)

    with open(args.output, "w") as f:
        json.dump(all_results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
