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
- `cage_references_ensembl.tsv`: reference manifest with explicit Ensembl FASTA/GTF URLs;
- `species_manifest_index.tsv` and `species/<species>/`: one sample/reference
  manifest pair per covered species for independent runs;
- `ensembl_coverage_report.tsv`: species/run coverage summary;
- `manifest_build_summary.json`: census totals and generation evidence.

The committed reference manifest uses explicit Ensembl FASTA/GTF URLs and sets
`chrom_sizes` to `-`; chromosome sizes are derived from the downloaded FASTA.
For a private or pinned mirror, replace the URLs and optionally add
`chrom_sizes_url` and `*_checksum` columns. The pipeline downloads files into
`--reference_cache`, verifies any supplied checksums, and reuses matching cached
files.

The 2026-11-11 snapshot contains 18,658 ENA runs across 101 reported species:

- 18,196 runs across 61 species with an Ensembl assembly;
- 462 runs across 40 species/ambiguous taxa without Ensembl assembly coverage.

The covered production run uses:

```bash
nextflow run /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/pipelines/cage/main.nf \
  -profile slurm \
  -c /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/config/singularity.config \
  --sample_manifest /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/pipelines/cage/manifests/20261111/cage_samples_ensembl.tsv \
  --reference_manifest /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/pipelines/cage/manifests/20261111/cage_references_ensembl.tsv \
  --reference_cache /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/reference-cache \
  --outdir /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/results/20261111 \
  -work-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/work/20261111 \
  -resume
```

Review `ensembl_coverage_report.tsv` and `cage_samples_uncovered.tsv` before
deciding whether to add non-Ensembl references for excluded taxa.

The combined manifest is the canonical input for an all-species run. Use the
per-species pair when scheduling or rerunning one species independently.
