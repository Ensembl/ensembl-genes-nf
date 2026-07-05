# Ribo-seq Pipeline

The Ribo-seq pipeline processes ribosome profiling data from sample acquisition through read preparation, QC, alignment, RiboMetric/RiboWaltz analysis, and genome track generation.

## Entry Point

```bash
nextflow run pipelines/riboseq/main.nf \
  -c pipelines/riboseq/nextflow.config \
  -profile slurm \
  --sample_sheet samples.csv \
  --outdir results/riboseq \
  --star_index /path/to/star_index \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --chrom_sizes_file genome.chrom.sizes
```

`--fasta` is checked explicitly in `main.nf` because it is required for RiboWaltz P-site analysis.

## Sample Sheet

The sample sheet is a CSV with at least a `Run` column. `study_accession` is optional and defaults to `unknown` when absent.

```text
Run,study_accession
SRR000001,PRJEB00000
```

When `--fetch true`, runs that do not already have collapsed reads are downloaded and processed. When `--fetch false`, the pipeline expects collapsed reads to be available through `--collapsed_read_path`.

## Workflow Stages

| Stage | Subworkflow | Key outputs |
| --- | --- | --- |
| Data acquisition | `DATA_ACQUISITION` | Collapsed FASTA reads, FASTQ/QC reports when fetching. |
| Quality control | `QUALITY_CONTROL` | getRPF cleanliness reports. |
| Alignment | `ALIGNMENT` | Genome BAMs, transcriptome BAMs, STAR logs. |
| Analysis | `ANALYSIS` | RiboMetric reports, RiboWaltz offsets and QC tables. |
| Post-processing | `POST_PROCESSING` | Filtered BAMs, bedGraphs, BigWigs, optional unique-read outputs. |

## Main Output Directories

| Directory | Contents |
| --- | --- |
| `fastq/`, `fastqc/`, `fastp/` | Downloaded reads and read-preparation QC. |
| `collapsed_fa/` | Collapsed reads used downstream. |
| `star_align/`, `samtools_sort/` | Alignment and sorted BAM files. |
| `RiboMetric/`, `ribowaltz/` | QC reports, offset tables, and plots. |
| `filtered_bams/`, `bedgraphs/`, `bigwigs/` | Track-ready outputs. |
| `unique_reads/` | Optional unique-read tracking outputs. |

## Parameters

```{toctree}
:maxdepth: 1

../generated/riboseq-parameters
```
