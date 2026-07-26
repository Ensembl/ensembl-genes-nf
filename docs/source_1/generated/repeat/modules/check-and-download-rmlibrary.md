# checks

regarding copyright ownership.
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Process Details

| Property | Value |
|----------|-------|
| Process | `checks` |
| Label | `'default'` |
| Tag | `$meta.gca` |
| Publish directory | `"${params.outdir}/${meta.gca}/rm_library", mode: 'copy'` |

## Inputs

```text
tuple val(url), val(meta)
```

## Outputs

```text
tuple val(meta), path("*.fa"), emit: repeatmodeler_library_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Download external resources
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/check_and_download_rmlibrary.nf`
