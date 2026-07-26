# GENERATE_HTML_REPORT

## Process Details

| Property | Value |
|----------|-------|
| Process | `GENERATE_HTML_REPORT` |
| Label | `'process_high'` |
| Publish directory | `"${params.outdir}", mode: 'copy'` |
| container | `"oras://community.wave.seqera.io/library/pip_jinja2_pandas:fee727bf7c211ccf"` |

## Inputs

### Nextflow interface

```nextflow
path(consensus_results)
val(ucsc_session_url)
```

## Outputs

### Nextflow interface

```nextflow
path("translon_consensus_report.html"), emit: html_report
```

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/translon-consensus/modules/generate_report.nf`
