# COMBINE_OUTPUTS

## Process Details

| Property | Value |
|----------|-------|
| Process | `COMBINE_OUTPUTS` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/combined",` |
| container | `"https://depot.galaxyproject.org/singularity/ubuntu:20.04"` |
| errorStrategy | `'ignore'` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(input_files)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("${meta.id}_combined.txt"),   emit: combined
path "versions.yml",                                emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/example/modules/combine_outputs.nf`
