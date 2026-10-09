# ENA Alignment Evidence Submission

## Purpose

The ENA pipeline publishes processed RNA-seq alignments from an Ensembl
genebuild as annotation evidence. It translates the layout used by Ensembl
into the objects expected by ENA, while keeping the submission reproducible
and safe to test before production use.

The important design rule is:

```text
one assembly + partial release
  └── one ENA Project/Study
        └── one ENA ANALYSIS per BAM/CRAM alignment
```

ENA does not accept several alignment files in one `ANALYSIS`. The pipeline
therefore starts with one annotation-level manifest, expands it into one file
level record per alignment, and generates one analysis XML per alignment.

## Data flow

1. Build or provide an annotation manifest and its file manifest.
2. Prepare the INSDC reference and validate alignment headers.
3. Reheader BAMs or convert them to CRAM when requested.
4. Calculate checksums and indexes for the files that will actually be sent.
5. Generate project and analysis XML, upload files to Webin, and poll for
   accessions.
6. Publish XML, accession, execution-report, and software-version outputs.

The pipeline deliberately separates metadata generation from alignment
preparation. This makes retries predictable and means that a failed upload can
be repeated without rebuilding unrelated inputs.

## Start here

For the short command reference and parameter table, see the pipeline
[README](https://github.com/Ensembl/ensembl-genes-nf/blob/dev/pipelines/ena_submit/README.md)
in the repository.

For a complete first-time run, including credentials, TEST/PROD promotion,
input discovery, and troubleshooting, see the
[submission walkthrough](https://github.com/Ensembl/ensembl-genes-nf/blob/dev/pipelines/ena_submit/WALKTHROUGH.md).

For the design decisions and remaining policy questions, see the
[submission strategy](https://github.com/Ensembl/ensembl-genes-nf/blob/dev/pipelines/ena_submit/SUBMISSION_STRATEGY.md).

## Safety rules

- Run against ENA `test` before `prod`.
- Store `ENA_WEBIN_PASSWORD` in Nextflow secrets; do not pass it as a CLI
  parameter.
- Use a separate output directory for each submission.
- Inspect generated XML and `accessions.tsv` before treating a submission as
  complete.
