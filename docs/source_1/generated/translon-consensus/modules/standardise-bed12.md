# STANDARDISE_BED12

## Process Details

| Property | Value |
|----------|-------|
| Process | `STANDARDISE_BED12` |
| Label | `'process_high'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/standardised_bed12s/${tool}", mode: 'copy'` |
| container | `"oras://community.wave.seqera.io/library/pybedtools_pysam_pip_biopython:1dbd8151223e518c"` |

## Inputs

```text
tuple val(meta), val(tool), path(bed_file)
path genome_fasta
path genome_fasta_fai
```

## Outputs

```text
tuple val(meta), val(tool), path("*.valid.bed12")
```

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/translon-consensus/modules/standardise_bed12.nf`
