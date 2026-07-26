# FASTP

## Process Details

| Property | Value |
|----------|-------|
| Process | `FASTP` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/fastp", mode: 'copy', pattern: '*_clipped.fastq.gz'` |
| conda | `"bioconda::fastp=0.23.4"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

```text
tuple val(meta), path(reads), path(adapter_fasta)
```

## Outputs

```text
tuple val(meta), path("*_clipped_provided.fastq.gz"), emit: trimmed_fastq_provided
tuple val(meta), path("*_clipped_final.fastq.gz"), emit: trimmed_fastq
tuple val(meta), path("*_provided_fastp.json"), emit: json_provided
tuple val(meta), path("*_auto_fastp.json"), emit: json_auto
tuple val(meta), path("*_provided_fastp.html"), emit: html_provided
tuple val(meta), path("*_auto_fastp.html"), emit: html_auto
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/fastp.nf`
