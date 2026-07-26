# TOOL_A

## Process Details

| Property | Value |
|----------|-------|
| Process | `TOOL_A` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/tool_a",` |
| container | `"https://depot.galaxyproject.org/singularity/ubuntu:20.04"` |
| errorStrategy | `'ignore'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(input_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("${meta.id}_A.txt"),          emit: results
path "versions.yml",                                emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/example/modules/tool_a.nf`
