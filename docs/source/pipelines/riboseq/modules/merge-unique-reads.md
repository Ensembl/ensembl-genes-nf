# MERGE_UNIQUE_READS

## Process Details

| Property | Value |
|----------|-------|
| Process | `MERGE_UNIQUE_READS` |
| Label | `'process_high'` |
| Tag | `${mode}` |
| Publish directory | `"${params.outdir}/unique_reads", mode: 'copy', pattern: "processing_summary*"` |
| conda | `"conda-forge::python=3.10 conda-forge::polars=0.20.0"` |
| container | `'community.wave.seqera.io/library/pip_polars:50bb6fae7997c472'` |

## Inputs

### Nextflow interface

```nextflow
path collapsed_fastas       // List of all collapsed FASTA files from current run
path previous_index         // Optional: previous unique_reads.fa, count_matrix.parquet, read_mapping.json
val mode                    // 'per_run' or 'progressive'
```

## Outputs

### Nextflow interface

```nextflow
path "unique_reads*.fa", emit: unique_reads_fasta
path "count_matrix*.parquet", emit: count_matrix
path "read_mapping*.json", emit: read_mapping
path "processing_summary*.json", emit: summary
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/merge_unique_reads.nf`
