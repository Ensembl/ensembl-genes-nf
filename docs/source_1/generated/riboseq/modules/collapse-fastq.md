# COLLAPSE_FASTQ

## Process Details

| Property | Value |
|----------|-------|
| Process | `COLLAPSE_FASTQ` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/collapsed_fa", mode: 'copy', pattern: "*collapsed.fa"` |
| conda | `"conda-forge::python=3.10 pip::riboseq-dp-tools=0.1.10"` |
| container | `"ghcr.io/lapti-ucc/riboseqorg-nf-rdp-tools:latest"` |

## Inputs

```text
tuple val(meta), path(fastq)
```

## Outputs

```text
tuple val(meta), path("*collapsed.fa"), emit: collapsed_fasta
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/collapse_fastq.nf`
