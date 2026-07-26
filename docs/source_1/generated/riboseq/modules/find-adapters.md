# FIND_ADAPTERS

## Process Details

| Property | Value |
|----------|-------|
| Process | `FIND_ADAPTERS` |
| Label | `'process_low'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/adapter_reports", mode: 'copy', pattern: "*_adapter_report.fa"` |
| conda | `"conda-forge::python=3.11 conda-forge::biopython conda-forge::pandas"` |
| container | `"ghcr.io/lapti-ucc/riboseqorg-nf-python-pandas-sqlite:latest"` |

## Inputs

```text
tuple val(meta), path(raw_fastq)
tuple val(meta2), path(fastqc_data)
```

## Outputs

```text
tuple val(meta), path("*_adapter_report.fa"), emit: adapter_report
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/find_adapters.nf`
