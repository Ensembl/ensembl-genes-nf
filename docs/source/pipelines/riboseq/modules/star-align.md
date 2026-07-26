# STAR_ALIGN

## Process Details

| Property | Value |
|----------|-------|
| Process | `STAR_ALIGN` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `path: "${params.outdir}/star_align", mode: 'copy', saveAs: {` |
| conda | `"bioconda::star=2.7.11b"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(reads)
path index
path gtf
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
tuple val(meta), path("*.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
tuple val(meta), path("*.Log.final.out"), emit: log
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/star_align.nf`
