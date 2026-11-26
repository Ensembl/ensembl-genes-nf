# UPLOAD_INTO_FTP

## Process Details

| Property | Value |
|----------|-------|
| Process | `UPLOAD_INTO_FTP` |
| Label | `'ensembl_ftp'` |
| Tag | `${meta.gca}:upload_into_ftp` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(files, arity: '1..*')
val ftp_base
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path(files), emit: library_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/upload_into_ftp.nf`
