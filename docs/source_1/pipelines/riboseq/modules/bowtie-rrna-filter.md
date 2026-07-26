# BOWTIE_RRNA_FILTER

## Process Details

| Property | Value |
|----------|-------|
| Process | `BOWTIE_RRNA_FILTER` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/rrna_filter", mode: 'copy', pattern: "*_rrna_filter.log"` |
| conda | `"bioconda::bowtie=1.3.1 bioconda::samtools=1.19"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(reads)
path index  // Bowtie index directory or files
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*_no_rrna.fastq.gz"), emit: filtered_fastq
tuple val(meta), path("*_rrna_filter.log"), emit: log
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/bowtie_rrna_filter.nf`
