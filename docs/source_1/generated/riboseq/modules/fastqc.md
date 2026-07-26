# FASTQC

## Process Details

| Property | Value |
|----------|-------|
| Process | `FASTQC` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/fastqc", mode: 'copy', saveAs: { filename ->` |
| conda | `"bioconda::fastqc=0.12.1"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

```text
tuple val(meta), path(fastq)
path adapter_list
```

## Outputs

```text
tuple val(meta), path("*_fastqc.html"), emit: html
tuple val(meta), path("*_fastqc.zip"), emit: zip
tuple val(meta), path("*/*_fastqc_data.txt"), emit: txt
path "versions.yml", emit: versions
```

## Implementation Summary

- Rename output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/fastqc.nf`
