# fetches

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
| Process | `fetches` |
| Label | `'fetch_file'` |
| Tag | `$meta.gca:repeatmodel` |
| Publish directory | `"${params.outdir}/${meta.gca}/library/", mode: 'copy'` |

## Inputs

```text
val(meta)
```

## Outputs

```text
tuple val(meta), path("${meta.gca}.repeatmodeler.fa"), emit: rep_library_file_output
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Run RepeatModeler
- Generate software version report
- Download external resources

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/fetch_repeat_model.nf`
