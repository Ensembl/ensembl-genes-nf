# FASTQ_DL

## Process Details

| Property | Value |
|----------|-------|
| Process | `FASTQ_DL` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/fastq", mode: 'copy', pattern: '*.fastq.gz'` |
| conda | `"bioconda::fastq-dl=2.0.1"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), val(run), path(needs_processing)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.fastq.gz"), emit: fastq
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/fastq_dl.nf`
