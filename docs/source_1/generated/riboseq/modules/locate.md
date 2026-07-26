# LOCATE

## Process Details

| Property | Value |
|----------|-------|
| Process | `LOCATE` |
| Label | `'process_light'` |
| Tag | `${meta.id}` |

## Inputs

```text
tuple val(meta), val(run)
```

## Outputs

```text
tuple val(meta), path("${run}.collapsed.fa.gz"), emit: collapsed_reads, optional: true
tuple val(meta), val(run), path("${run}_needs_processing"), emit: needs_processing, optional: true
```

## Implementation Summary

- Create symbolic links

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/locate.nf`
