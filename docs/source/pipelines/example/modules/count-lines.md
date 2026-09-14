# COUNT_LINES

## Process Details

| Property | Value |
|----------|-------|
| Process | `COUNT_LINES` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/line_counts",` |
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
tuple val(meta), path("${meta.id}_line_count.txt"), emit: counts
path "versions.yml",                                emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/example/modules/count_lines.nf`
