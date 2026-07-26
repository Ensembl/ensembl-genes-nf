# MERGE_BEDGRAPHS

## Process Details

| Property | Value |
|----------|-------|
| Process | `MERGE_BEDGRAPHS` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/bedgraphs", mode: 'copy', pattern: "*.merged.*.sorted.bedgraph"` |
| conda | `"bioconda::bedtools=2.31.1"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

```text
tuple val(meta), path(bedgraphs)
```

## Outputs

```text
tuple val(meta), path("*.merged.*.sorted.bedgraph"), emit: bedgraph
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/merge_bedgraphs.nf`
