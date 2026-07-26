# TOOL_D

## Process Details

| Property | Value |
|----------|-------|
| Process | `TOOL_D` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/tool_d",` |
| container | `"https://depot.galaxyproject.org/singularity/ubuntu:20.04"` |
| errorStrategy | `'ignore'` |

## Inputs

```text
tuple val(meta), path(input_file)
```

## Outputs

```text
tuple val(meta), path("${meta.id}_D.txt"),          emit: results
path "versions.yml",                                emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/example/modules/tool_d.nf`
