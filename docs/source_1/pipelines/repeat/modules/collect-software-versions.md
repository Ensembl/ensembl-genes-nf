# to

regarding copyright ownership.
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
distributed under the License is distributed on an "AS IS" BASIS,
See the License for the specific language governing permissions and
limitations under the License.

## Process Details

| Property | Value |
|----------|-------|
| Process | `to` |
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
