# uploads

regarding copyright ownership.
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
distributed under the License is distributed on an "AS IS" BASIS,
See the License for the specific language governing permissions and
limitations under the License.

## Process Details

| Property | Value |
|----------|-------|
| Process | `uploads` |
| Label | `'ensembl_ftp'` |
| Tag | `$meta.gca:upload_repeats_into_ftp` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta),path(fasta_file),path(stk_file),path(log_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta),path(fasta_file), emit: library_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/upload_repeats_into_ftp.nf`
