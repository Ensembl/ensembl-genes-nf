# RIBOWALTZ

## Process Details

| Property | Value |
|----------|-------|
| Process | `RIBOWALTZ` |
| Label | `'process_medium'` |
| Tag | `${meta.id}` |
| Publish directory | `"${params.outdir}/ribowaltz", mode: 'copy', pattern: "*.{tsv.gz,pdf}", saveAs: { filename -> filename.endsWith('.pdf') ? "offset_plot/${filename}" : filename }` |
| conda | `"bioconda::bioconductor-ribowaltz=2.0"` |
| container | `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?` |

## Inputs

```text
tuple val(meta), path(transcriptome_bam), path(transcriptome_bam_index)
tuple val(meta2), path(gtf)
tuple val(meta3), path(fasta)
```

## Outputs

```text
tuple val(meta), path("*.cds_coverage_psite.tsv.gz"), optional: true, emit: cds_coverage
tuple val(meta), path("offset_plot"), optional: true, emit: offset_plots
tuple val(meta), path("*.psite_offset.tsv.gz"), optional: true, emit: psite_offsets
tuple val(meta), path("*.psite.tsv.gz"), optional: true, emit: psite_table
tuple val(meta), path("*nt_coverage_psite.tsv.gz"), optional: true, emit: nt_coverage
tuple val(meta), path("*.codon_coverage_rpf.tsv.gz"), optional: true, emit: codon_rpf
tuple val(meta), path("*.codon_coverage_psite.tsv.gz"), optional: true, emit: codon_psite
tuple val(meta), path("ribowaltz_qc/*.pdf"), optional: true, emit: qc_plots
tuple val(meta), path("*.best_offset.txt"), optional: true, emit: best_offset
path "versions.yml", emit: versions
```

## Implementation Summary

- Generate software version report

## Source

`/Users/ftricomi/Documents/ensembl-genes-nf/pipelines/riboseq/modules/ribowaltz.nf`
