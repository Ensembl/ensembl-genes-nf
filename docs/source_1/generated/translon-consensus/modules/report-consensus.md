# to

## Process Details

| Property | Value |
|----------|-------|
| Process | `to` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/consensus_reports/${meta.id}", mode: 'copy'` |
| container | `"oras://community.wave.seqera.io/library/pip_duckdb_pandas:57c3741f55d53490"` |

## Inputs

```text
tuple val(meta), path(samplesheet), path(bed_files)
val ucsc_session_url
path gencode_gtf, stageAs: 'gencode.gtf*'
```

## Outputs

```text
tuple val(meta), path("*.tsv"),           emit: results
```

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/translon-consensus/modules/report_consensus.nf`
