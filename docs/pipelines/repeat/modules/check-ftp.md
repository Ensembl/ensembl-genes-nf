# CHECK_FTP

This process fetches the RepeatModeler library file for a given GCA accession.
If the file is available at the specified URL, it downloads the file
and saves it with the name "<GCA>.repeatmodeler.fa".
If the file is not available, it outputs a message indicating that
the download was skipped and creates an empty file with the same name.
The fetched library file is saved in the "library" directory under
the output directory for the given GCA accession.

## Process Details

| Property | Value |
|----------|-------|
| Process | `CHECK_FTP` |
| Label | `'fetch_file'` |
| Tag | `${meta.gca}:check_${file_type}` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), val(url), val(file_type)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), val(url), val(file_type), path("status.txt"),
emit: status_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report
- Download external resources

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/check_ftp.nf`
