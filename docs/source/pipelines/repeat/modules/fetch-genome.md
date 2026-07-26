# FETCH_GENOME

regarding copyright ownership.
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
distributed under the License is distributed on an "AS IS" BASIS,
See the License for the specific language governing permissions and
limitations under the License.

## Process Details

| Property | Value |
|----------|-------|
| Process | `FETCH_GENOME` |
| Label | `'python'` |
| Tag | `${meta.gca}:genome` |
| storeDir | `"${params.cacheDir}/${meta.gca}/ncbi_dataset"` |
| maxForks | `10` |

## Inputs

### Nextflow interface

```nextflow
val meta
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.fa"), emit: genome_file_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Copy output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/fetch_genome.nf`
