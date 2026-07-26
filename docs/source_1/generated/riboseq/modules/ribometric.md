# RIBOMETRIC

## Process Details

| Property | Value |
|----------|-------|
| Process | `RIBOMETRIC` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/RiboMetric", mode: 'copy', pattern: "*RiboMetric.{html,json,csv}"` |
| conda | `"conda-forge::python=3.10 conda-forge::biopython bioconda::pysam"` |
| container | `"ghcr.io/lapti-ucc/riboseqorg-nf-ribometric:latest"` |
| errorStrategy | `'ignore'` |

## Inputs

```text
tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
path ribometric_annotation
```

## Outputs

```text
tuple val(meta), path("*RiboMetric.html"), emit: html
tuple val(meta), path("*RiboMetric.json"), emit: json
tuple val(meta), path("*RiboMetric.csv"), emit: csv
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/ribometric.nf`
