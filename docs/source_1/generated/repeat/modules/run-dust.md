# runs

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
| Process | `runs` |
| Label | `python` |
| Tag | `${meta.gca}:genome` |
| Publish directory | `"${params.outdir}/${meta.gca}/dust/", pattern: "**/*.gtf", mode: "copy"` |

## Inputs

```text
val(meta)
```

## Outputs

```text
tuple val(meta), path("*.gtf"), emit: dust_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- run dust
- Rename output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/run_dust.nf`
