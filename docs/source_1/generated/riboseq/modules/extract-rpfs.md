# EXTRACT_RPFS

## Process Details

| Property | Value |
|----------|-------|
| Process | `EXTRACT_RPFS` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/getRPF/extract", mode: 'copy', pattern: "*.{seqspec.yaml,extraction_report.json,report.html}"` |
| conda | `"conda-forge::python=3.10 conda-forge::biopython"` |
| container | `"ghcr.io/jackcurragh/get-rpf:main"` |

## Inputs

```text
tuple val(meta), path(input_file)
path star_index
```

## Outputs

```text
tuple val(meta), path("*_rpfs.fastq"), emit: rpfs
tuple val(meta), path("*.seqspec.yaml"), emit: seqspec
tuple val(meta), path("*.extraction_report.json"), emit: report
tuple val(meta), path("*.report.html"), emit: html_report
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/extract_rpfs.nf`
