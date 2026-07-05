# Translon Consensus Pipeline

The translon consensus pipeline compares ORF predictions from multiple tools and produces per-sample consensus tables plus an HTML summary report.

## Entry Point

```bash
nextflow run pipelines/translon-consensus/main.nf \
  -c pipelines/translon-consensus/nextflow.config \
  --bed_results_dir results_bed12 \
  --gencode_gtf gencode.v45.annotation.gtf.gz \
  --gencode_fasta genome.fa \
  --gencode_fasta_fai genome.fa.fai \
  --outdir results/translon-consensus
```

The pipeline expects BED or BED12 files grouped by tool:

```text
results_bed12/
├── RiboTIE/
│   └── sample_a.bed12
├── PRICE/
│   └── sample_a.bed12
└── ORFQuant/
    └── sample_a.bed12
```

Samples with only one tool are logged and skipped because consensus analysis needs at least two callers.

## Preparing BED12 Input

Tool-specific converters live in `pipelines/translon-consensus/scripts/`. The recommended wrapper auto-detects common formats:

```bash
python pipelines/translon-consensus/scripts/convert_to_bed12.py input.bed results_bed12/RiboTIE/
```

Use the individual converters when auto-detection is not enough:

- `ribotie_to_bed12.py`
- `price_to_bed12.py`
- `iribo_to_bed12.py`
- `orfquant_to_bed12.py`

## Workflow Stages

| Stage | Module | Purpose |
| --- | --- | --- |
| Standardise BED12 | `STANDARDISE_BED12` | Validate and normalise BED/BED12 inputs. |
| Rename by tool | `RENAME_BED` | Prefix files with the calling tool name. |
| Create samplesheet | `CREATE_SAMPLESHEET` | Build per-sample inputs for consensus analysis. |
| Consensus report | `REPORT_CONSENSUS` | Produce per-sample consensus TSV files. |
| HTML report | `GENERATE_HTML_REPORT` | Build `translon_consensus_report.html`. |
| Optional FTP move | `MOVE_TO_FTP` | Copy the HTML report to `--ftp_dir` when supplied. |

## Outputs

| Directory or file | Contents |
| --- | --- |
| `standardised_bed12s/<tool>/` | Valid BED12 files grouped by tool. |
| `consensus_reports/<sample>/` | Consensus TSV outputs for multi-tool samples. |
| `translon_consensus_report.html` | Interactive HTML summary report. |

## Parameters

```{toctree}
:maxdepth: 1

../generated/translon-consensus-parameters
```
