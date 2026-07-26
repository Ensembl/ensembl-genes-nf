# SAMTOOLS_SORT

## Process Details

| Property | Value |
|----------|-------|
| Process | `SAMTOOLS_SORT` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/samtools_sort", mode: 'copy', pattern: "*.sorted.bam"` |
| conda | `"bioconda::samtools=1.20"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(bam)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.sorted.bam"), emit: bam
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/samtools_sort.nf`
