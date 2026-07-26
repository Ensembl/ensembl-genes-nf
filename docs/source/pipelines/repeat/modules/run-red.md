# RUN_RED

regarding copyright ownership.
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
distributed under the License is distributed on an "AS IS" BASIS,
See the License for the specific language governing permissions and
limitations under the License.

## Process Details

| Property | Value |
|----------|-------|
| Process | `RUN_RED` |
| Label | `python` |
| Tag | `${meta.gca}:genome` |
| Publish directory | `"${params.outdir}/${meta.gca}/red/", pattern: "**/*.gtf", mode: "copy"` |

## Inputs

### Nextflow interface

```nextflow
val(meta)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.gtf"), emit: red_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Execute Run Red
- Rename output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/run_red.nf`
