# BAM_TO_BED

## Process Details

| Property | Value |
|----------|-------|
| Process | `BAM_TO_BED` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/bedgraphs", mode: 'copy', pattern: "*.sorted.bedgraph"` |
| conda | `"conda-forge::python=3.9 bioconda::pysam=0.23.3 bioconda::samtools=1.20 conda-forge::biopython conda-forge::numpy"` |
| container | `"ghcr.io/lapti-ucc/riboseqorg-nf-bam-to-bed:latest"` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(bam), path(bai), path(offsets)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.sorted.bedgraph"), emit: bedgraph
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/bam_to_bed.nf`
