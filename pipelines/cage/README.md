# Protocol-aware CAGE atlas pipeline

This pipeline processes a curated CAGE or related 5-prime assay manifest against
native, versioned species references. It does not search or download the full
archive inventory implicitly. Use the survey inventory to create a reviewed
sample manifest, then run this pipeline explicitly.

## Manifests

`sample_manifest` is tab-separated with required columns:

```text
sample_id species assembly assay_class protocol layout accession fastq_1 fastq_2
```

`assay_class` is one of `cage`, `rampage`, `pro-cap`, `net-cage`, or `tss-seq`.
`layout` is `SE` or `PE`. Supply either an archive `accession` or local
`fastq_1`; paired-end rows also require `fastq_2` when local files are used.

`reference_manifest` is tab-separated with:

```text
assembly species fasta gtf chrom_sizes fasta_url gtf_url chrom_sizes_url taxon_id source
```

The `fasta`, `gtf`, and `chrom_sizes` columns may contain local files, or the
manifest may use `fasta_url`, `gtf_url`, and `chrom_sizes_url`. Optional
`fasta_checksum`, `gtf_checksum`, and `chrom_sizes_checksum` columns accept
SHA256 values. Downloads are stored in `params.reference_cache` and reused
when the source and checksum are unchanged. The sample assembly must match
one row in the reference manifest. During `prepare_references`, rRNA features
are extracted from the supplied GFF/GTF and added as decoy contigs to the STAR
reference. Reads aligned to those decoys are reported through the reference
outputs and excluded from the CTSS layer.

An Ensembl-first production snapshot can be generated with
`bin/build_ensembl_manifests.py`; see `manifests/README.md`. It produces a
covered sample manifest plus an explicit coverage report for survey species
without an Ensembl assembly.

## Run

```bash
nextflow run /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/pipelines/cage/main.nf \
  -profile test \
  -c /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/ensembl-genes-nf/config/singularity.config \
  --reference_cache /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/reference-cache \
  --outdir /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/results/test \
  -work-dir /hps/nobackup/flicek/ensembl/genebuild/jackt/cage/work/test \
  -resume

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

All production tools run in Singularity containers, including the `test`
profile. The test profile only reduces resources and uses the local executor;
it does not disable containers.

The pipeline preserves `meta` fields through acquisition, trimming, alignment,
CTSS extraction, and atlas construction. CTSS BED/bedGraph files are immutable
published intermediates. Related assay classes are labelled and are not pooled
silently.

## Outputs

FASTQs are staged internally for processing and are not published by default.
Published outputs are grouped under `alignments/`, `ctss/`, and
`atlas/<species>/<assembly>/`. The atlas directory contains sample metadata,
CTSS/TSS/promoter matrices, TSS clusters, bidirectional enhancer candidates,
GTF-overlap annotations, and QC metrics. The cross-species comparison is opt-in with a JSON params file
containing `"run_cross_species": true` and an `orthology_manifest` path.

The test profile is a wiring check. Production analysis requires validated
references, protocol-specific parameters, real container execution, and
biological QC.
