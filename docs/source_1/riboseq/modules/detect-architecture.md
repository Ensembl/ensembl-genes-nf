# DETECT_ARCHITECTURE

## Process Details

| Property | Value |
|----------|-------|
| Process | `DETECT_ARCHITECTURE` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/getRPF/detect", mode: 'copy', pattern: "*.{seqspec.yaml,adapters.fa,extraction_report.json,fastq}"` |
| conda | `"conda-forge::python=3.10 conda-forge::biopython"` |
| container | `"ghcr.io/lapti-ucc/riboseqorg-nf-getrpf:latest"` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(input_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.seqspec.yaml"), emit: seqspec
tuple val(meta), path("*.adapters.fa"), emit: adapters
tuple val(meta), path("*.extraction_report.json"), emit: report
tuple val(meta), path("*_rpfs.fastq"), emit: rpfs, optional: true
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/detect_architecture.nf`
