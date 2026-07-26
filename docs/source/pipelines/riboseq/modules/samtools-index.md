# SAMTOOLS_INDEX

## Process Details

| Property | Value |
|----------|-------|
| Process | `SAMTOOLS_INDEX` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| conda | `"bioconda::samtools=1.20"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(sorted_bam)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("${sorted_bam}"), path("*.bai"), emit: bam_and_bai
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/samtools_index.nf`
