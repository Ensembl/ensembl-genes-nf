# RUN_REPEATMASKER

This process runs RepeatMasker on a given genome file using a
specified RepeatModeler library. It uses the RepeatMasker tool
to perform the analysis and generates GTF files containing the identified
repetitive regions. The output GTF files are saved in the "repeatmasker"
directory under the output directory for the given GCA accession.
The process also generates a versions.yml file containing the version
of RepeatMasker used.

## Process Details

| Property | Value |
|----------|-------|
| Process | `RUN_REPEATMASKER` |
| Label | `python` |
| Tag | `${meta.gca}:genome` |
| Publish directory | `"${params.outdir}/${meta.gca}/repeatmasker/", pattern: "repeatmasker_output/*.gtf", mode: "copy"` |

## Inputs

### Nextflow interface

```nextflow
tuple val(meta), val(library_file)
```

## Outputs

### Nextflow interface

```nextflow
tuple val(meta), path("*.gtf"), emit: repeatmasker_out
path "versions.yml", emit: versions_file
```

## Implementation Summary

- Execute Run Repeatmasker
- Rename output files
- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/repeat/modules/run_repeatmasker.nf`
