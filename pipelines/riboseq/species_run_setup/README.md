# RiboSeq Ensembl Beta species run setup

Generated from `/Users/jackt/FINAL_LIST.csv` and Ensembl Beta FTP (`https://ftp.ebi.ac.uk/pub/ensemblorganisms/`).

## Contents

- `sample_sheets/*_samplesheet.csv`: pipeline-ready CSVs with exactly `Run,study_accession`.
- `configs/*_beta_params.config`: URL-based configs that use Ensembl Beta FTP files directly.
- `species_manifest.tsv`: run counts, TaxIDs, inhibitor summaries, chosen Beta species directory, assembly, provider, geneset date, FASTA URL, and GTF URL.

## Run pattern on Codon

```bash
nextflow run /hps/software/users/ensembl/genebuild/jackt/ensembl-genes-nf/pipelines/riboseq/main.nf \
  -c /homes/jackt/riboseq_species_runs/configs/<slug>_beta_params.config \
  -profile slurm
```

The configs intentionally use `download_method = "url"` and direct Beta FTP URLs. They do not use `gget` or old Ensembl release-number FTP paths.

## Default reference rules

- Use canonical/reference assembly where the species has a clear standard assembly, for example human GRCh38.p14, mouse GRCm39, and zebrafish GRCz11.
- Otherwise use the newest complete Beta `genes.gtf.gz` record discoverable for the selected species/strain directory.
- Use `genes.gtf.gz`, not `genes-including_alt.gtf.gz`, for clean cross-species defaults and to avoid alt-locus duplicate annotation in ordinary Ribo-seq runs.
- Species without a confident Beta mapping remain sample-sheet-only in the manifest.
