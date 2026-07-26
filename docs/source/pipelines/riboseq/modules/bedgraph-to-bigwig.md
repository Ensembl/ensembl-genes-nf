# BEDGRAPH_TO_BIGWIG

## Process Details

| Property | Value |
|----------|-------|
| Process | `BEDGRAPH_TO_BIGWIG` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/bigwigs", mode: 'copy', pattern: "*.bw"` |
| conda | `"bioconda::ucsc-bedgraphtobigwig=469"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), path(bedgraph)
path chrom_sizes
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.bw"), emit: bigwig
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/bedgraph_to_bigwig.nf`
