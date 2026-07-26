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
| Label | `'repeatmodeler'` |
| Tag | `$meta.gca:run_repeatmodeler` |
| Publish directory | `"${params.outdir}/${meta.gca}/library", mode: 'copy'` |

## Inputs

```text
val(meta)
```

## Outputs

```text
tuple val(meta), path("*repeatmodeler.fa"), path("*families.stk.gz"), path("*rmod.log"), emit: repeatmodeler_library_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Run RepeatModeler
- Rename output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/generate_repeatmodeler_library.nf`
