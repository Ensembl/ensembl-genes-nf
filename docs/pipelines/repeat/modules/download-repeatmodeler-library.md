# DOWNLOAD_REPEATMODELER_LIBRARY

This process checks if a RepeatModeler library file exists at the given URL.
If the file exists, it downloads the file and saves it with the name "<GCA>.
repeatmodeler.fa". If the file does not exist, it outputs an error message
and exits with a non-zero status.

## Process Details

| Property | Value |
|----------|-------|
| Process | `DOWNLOAD_REPEATMODELER_LIBRARY` |
| Label | `'default'` |
| Tag | `$meta.gca` |
| Publish directory | `"${params.outdir}/${meta.gca}/rm_library", mode: 'copy'` |

## Inputs

| Name | Description |
| ---- | ----------- |
| `meta` | Genome metadata map |
| `url` | RepeatModeler library URL |

### Nextflow interface

```nextflow
tuple val(url), val(meta)
```

## Outputs

| Name | Description |
| ---- | ----------- |
| `meta` | Genome metadata map |
| `fasta_file` | RepeatModeler library file |

### Nextflow interface

```nextflow
tuple val(meta), path("*.fa"), emit: repeatmodeler_library_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report
- Download external resources

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/download_repeatmodeler_library.nf`
