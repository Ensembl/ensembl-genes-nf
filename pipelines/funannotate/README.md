# Funannotate BAM-backed prediction

This pipeline runs one `funannotate predict` task per assembly/cohort manifest
row. The soft-masked genome is read from the working filesystem. CRAMs and
optional CRAI files are acquired from FTP, converted to BAM with samtools using
the local genome as reference, indexed when needed, and passed to funannotate.
Optional transcript, protein, PASA, and pre-trained parameter files are also supported.

The genome must be soft-masked. CRAM URLs, CRAI URLs, and checksum lists are
semicolon-separated and aligned by position. CRAI URLs may be omitted; the
pipeline creates BAM BAI files after conversion.

Optional FTP evidence columns include `stringtie_ftp_url`,
`protein_alignments_ftp_url`, and `trnascan_ftp_url`; these map to the
corresponding Funannotate options.

The minimum manifest is:

```text
cohort_id	assembly_accession	assembly_release	species	isolate	genome_path	cram_ftp_urls	crai_ftp_urls	cram_sha256s
```

`cram_ftp_urls`, `crai_ftp_urls`, and `cram_sha256s` are semicolon-separated
lists in matching order. `genome_path` is resolved on the filesystem where
Nextflow is launched and must be visible to the compute nodes.

Example:

```bash
nextflow run pipelines/funannotate/main.nf -profile slurm \
  --input_manifest /path/to/funannotate.tsv \
  --stage_dir /hps/nobackup/flicek/ensembl/genebuild/jackt/funannotate-stage \
  --funannotate_container 'quay.io/biocontainers/funannotate@sha256:PINNED_DIGEST' \
  --funannotate_db /path/to/funannotate-db \
  --outdir /path/to/results \
  -resume
```

For a local stub wiring check:

```bash
nextflow run pipelines/funannotate/main.nf -stub-run -profile test \
  --input_manifest pipelines/funannotate/test/data/funannotate_manifest.tsv \
  --stage_dir /tmp/funannotate-stage \
  --outdir /tmp/funannotate-stub
```
