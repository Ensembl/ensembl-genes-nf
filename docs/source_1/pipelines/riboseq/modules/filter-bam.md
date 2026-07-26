# FILTER_BAM

## Process Details

| Property | Value |
|----------|-------|
| Process | `FILTER_BAM` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/filtered_bams", mode: 'copy', pattern: "*.bam"` |
| conda | `"bioconda::pysam=0.23.3"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(bam)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("${prefix}.unique_no_junction.bam"), emit: unique_no_junction
tuple val(meta), path("${prefix}.unique_with_junction.bam"), emit: unique_with_junction
tuple val(meta), path("${prefix}.multi_no_junction.bam"), emit: multi_no_junction
tuple val(meta), path("${prefix}.multi_with_junction.bam"), emit: multi_with_junction
tuple val(meta), path("${prefix}*.bam"), emit: all_bams
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/filter_bam.nf`
