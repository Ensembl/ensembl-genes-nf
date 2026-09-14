# COLLECT_SOFTWARE_VERSIONS

## Process Details

| Property | Value |
|----------|-------|
| Process | `COLLECT_SOFTWARE_VERSIONS` |
| Publish directory | `"${params.outdir}/pipeline_info", mode: 'copy'` |

## Inputs

### Nextflow interface

```nextflow
path 'versions_*.yml'
```

## Outputs

### Nextflow interface

```nextflow
path "software_versions.yml"
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/collect_software_versions.nf`
