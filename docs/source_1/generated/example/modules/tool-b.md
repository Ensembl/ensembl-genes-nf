# TOOL_B

## Process Details

| Property | Value |
|----------|-------|
| Process | `TOOL_B` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/tool_b",` |
| container | `"https://depot.galaxyproject.org/singularity/ubuntu:20.04"` |
| errorStrategy | `'ignore'` |

## Inputs

```text
tuple val(meta), path(input_file)
```

## Outputs

```text
tuple val(meta), path("${meta.id}_B.txt"),          emit: results
path "versions.yml",                                emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/example/modules/tool_b.nf`
