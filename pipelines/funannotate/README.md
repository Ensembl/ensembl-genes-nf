# Funannotate BAM-backed prediction

This pipeline runs one `funannotate predict` task per assembly/cohort manifest
row. Genome FASTA and coordinate-sorted RNA-seq BAMs are acquired from FTP,
checksummed when requested, copied to a persistent staging directory, indexed
when needed, and passed to funannotate. Optional transcript, protein, PASA, and
pre-trained parameter files are also supported.

The genome must be soft-masked. BAM URLs, index URLs, and checksum lists are
semicolon-separated and aligned by position. Index URLs may be omitted; the
pipeline will create BAI files with samtools.

Optional FTP evidence columns include `stringtie_ftp_url`,
`protein_alignments_ftp_url`, and `trnascan_ftp_url`; these map to the
corresponding Funannotate options.

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
