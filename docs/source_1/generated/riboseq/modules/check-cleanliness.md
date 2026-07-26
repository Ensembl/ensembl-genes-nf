# CHECK_CLEANLINESS

## Process Details

| Property | Value |
|----------|-------|
| Process | `CHECK_CLEANLINESS` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/getRPF/check", mode: 'copy', pattern: "*_{report,rpf_checks}.txt"` |
| conda | `"conda-forge::python=3.10 conda-forge::biopython"` |
| container | `"ghcr.io/jackcurragh/get-rpf:main"` |

## Inputs

```text
tuple val(meta), path(input_file)
val count_pattern
```

## Outputs

```text
tuple val(meta), path("*_report.txt"), emit: report
tuple val(meta), path("*rpf_checks.txt"), emit: rpf_checks
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/check_cleanliness.nf`
