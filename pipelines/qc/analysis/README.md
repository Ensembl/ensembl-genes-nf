# Pairwise benchmark analysis

The Nextflow workflow produces the reproducible comparison tables. This
directory contains lightweight post-workflow analysis tools so plotting and
report iteration does not require rerunning or copying the per-pair outputs.

## Make a compact bundle on the cluster

```bash
QC_ROOT=/hps/nobackup/flicek/ensembl/genebuild/jackt/long-read/qc_pairwise_20261101
python3 \
  "$QC_ROOT/ensembl-genes-nf/pipelines/qc/bin/extract_pairwise_benchmark.py" \
  --results-dir "$QC_ROOT/results/full/qc/pairwise_benchmark" \
  --output-dir "$QC_ROOT/analysis/full_bundle"
```

The extractor streams the large aggregate TSVs and writes compact summaries:

```text
pair_summary.tsv
tool_technology_summary.tsv
model_inventory_summary.tsv
structural_feature_summary.tsv
whole_model_class_counts.tsv
metric_status_summary.tsv
bundle_manifest.json
```

If `AGGREGATE_PAIRWISE_BENCHMARK` has not finished yet, use the same command
against the per-pair benchmark directory. The extractor will discover completed
`*.models.tsv`, `*.comparisons.tsv`, `*.whole_models.tsv`, and `*.summary.tsv`
files and mark the bundle as a partial run:

```bash
python3 \
  "$QC_ROOT/ensembl-genes-nf/pipelines/qc/bin/extract_pairwise_benchmark.py" \
  --results-dir "$QC_ROOT/results/full/qc/pairwise_benchmark" \
  --output-dir "$QC_ROOT/analysis/partial_bundle"
```

For an immediate first look that reads only the small per-pair summary files,
add `--summary-only`. This does not scan the large model or comparison tables:

```bash
python3 \
  "$QC_ROOT/ensembl-genes-nf/pipelines/qc/bin/extract_pairwise_benchmark.py" \
  --results-dir "$QC_ROOT/results/full/qc/pairwise_benchmark" \
  --output-dir "$QC_ROOT/analysis/summary_bundle" \
  --summary-only
```

No transcript-level tables are copied by default. To export detail for one or
more pairs for a focused investigation:

```bash
python3 \
  "$QC_ROOT/ensembl-genes-nf/pipelines/qc/bin/extract_pairwise_benchmark.py" \
  --results-dir "$QC_ROOT/results/full/qc/pairwise_benchmark" \
  --output-dir "$QC_ROOT/analysis/hg002a_iso_cDNA_bundle" \
  --pair-id hg002_HG002a_ISO_cDNA_bambu_vs_flair
```

Only the selected pairs' model, structural-comparison, and whole-model rows
are written under `details/`. This makes the bundle suitable for `rsync` back
to a local machine without moving the full result tree.

## Run the notebook on the cluster

```bash
export PAIRWISE_BUNDLE="$QC_ROOT/analysis/full_bundle"
jupyter lab --no-browser --ip=127.0.0.1
```

Open `pairwise_benchmark_analysis.ipynb`. It uses only the compact bundle and
writes plots to `$PAIRWISE_BUNDLE/plots/`. The notebook is deliberately
descriptive: it visualises model production, structural concordance, whole-
model classes, and metric status. It does not label pairwise agreement as
recall or precision.

## What still requires a truth set

The bundle can answer how tools differ in model counts, loci, intron chains,
ends, UTR structures, whole-model concordance, and evidence support. Recall,
precision, and F1 require a declared read-derived or independent reference
truth set and should be added as a separate analysis layer.
