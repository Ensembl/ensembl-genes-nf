# FETCH_GENOME

This process fetches the genome file for a given GCA accession.
If a genome file is provided as part of the metadata,
it uses that file instead of downloading it.
The fetched genome file is saved with the name "genome.fna".

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
