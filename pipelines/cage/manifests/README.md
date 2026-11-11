# Production manifest generation

The production manifests are generated from the current ENA structured CAGE
census and Ensembl's current species catalogue. The generator does not download
FASTQs or reference genomes.

```bash
python pipelines/cage/bin/build_ensembl_manifests.py \
  --output-dir pipelines/cage/manifests/20261111 \
  --reference-root /path/to/ensembl-cage-references
```

Outputs:

- `cage_samples_ensembl.tsv`: ENA runs whose species have an Ensembl assembly;
- `cage_samples_uncovered.tsv`: excluded runs without Ensembl assembly coverage;
- `cage_references_ensembl.tsv`: reference manifest with paths to be staged;
- `species_manifest_index.tsv` and `species/<species>/`: one sample/reference
  manifest pair per covered species for independent runs;
- `ensembl_coverage_report.tsv`: species/run coverage summary;
- `manifest_build_summary.json`: census totals and generation evidence.

The reference manifest is intentionally marked `needs_staging` until the FASTA,
GTF, and chromosome-size files are placed at the configured reference root.

The 2026-11-11 snapshot contains 18,658 ENA runs across 101 reported species:

- 18,196 runs across 61 species with an Ensembl assembly;
- 462 runs across 40 species/ambiguous taxa without Ensembl assembly coverage.

The covered production run uses:

```bash
nextflow run pipelines/cage/main.nf -profile slurm \
  --sample_manifest pipelines/cage/manifests/20261111/cage_samples_ensembl.tsv \
  --reference_manifest pipelines/cage/manifests/20261111/cage_references_ensembl.tsv \
  --outdir /path/to/cage-results -resume
```

Review `ensembl_coverage_report.tsv` and `cage_samples_uncovered.tsv` before
deciding whether to add non-Ensembl references for excluded taxa.

The combined manifest is the canonical input for an all-species run. Use the
per-species pair when scheduling or rerunning one species independently.
