# CREATE_SAMPLESHEET

## Process Details

| Property | Value |
|----------|-------|
| Process | `CREATE_SAMPLESHEET` |
| Tag | `${meta.id}` |

## Inputs

```text
tuple val(meta), path(bed_files)
```

## Outputs

```text
tuple val(meta), path("${meta.id}_samplesheet.tsv"), path(bed_files)
```

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/translon-consensus/modules/create_samplesheet.nf`
