# RENAME_BED

## Process Details

| Property | Value |
|----------|-------|
| Process | `RENAME_BED` |
| Label | `'process_light'` |
| Tag | `${meta.id} - ${tool}` |

## Inputs

```text
tuple val(meta), val(tool), path(bed_file)
```

## Outputs

```text
tuple val(meta), path("${tool}_${bed_file.name}")
```

## Implementation Summary

- Copy generated files

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/translon-consensus/modules/rename_bed.nf`
