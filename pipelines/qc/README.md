# QC Pipeline

This pipeline runs annotation quality-control metrics for one or more GFF3 files.

At the moment, the implemented QC step is AGAT statistics generation followed by conversion of the AGAT text report into a `genebuild` CSV using the `parse_agat.py` parser from the `ensembl-genes` repository.

## What It Does

For each row in an input CSV, the pipeline:

1. runs `agat_sp_statistics.pl` on the input GFF3 file
2. writes the AGAT text report
3. parses that report into a `*_genebuild.csv` metrics file

The workflow processes all listed samples in parallel.

## Requirements

- Nextflow with DSL2 enabled
- Singularity or Apptainer available on the execution host
- A local checkout of `ensembl-genes` containing:
  - `src/python/ensembl/genes/annotation-qc/parsers/parse_agat.py`
  - `src/python/ensembl/genes/annotation-qc/metrics/config/feature_levels.yaml`

## Inputs

### Sample Sheet

The required input is `--gff_csv`, a CSV file with a header and at least these columns:

```csv
sample,gff3
sample_1,/path/to/sample_1.gff3
sample_2,/path/to/sample_2.gff3
```

- `sample`: sample identifier used for logging and task tags
- `gff3`: absolute or relative path to the GFF3 file

## Parameters

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--gff_csv` | yes | none | CSV describing the samples and GFF3 paths |
| `--ensembl_genes_repo` | yes | none | Path to a local `ensembl-genes` checkout |
| `--feature_levels` | no | derived from `ensembl_genes_repo` | Path to `feature_levels.yaml` |
| `--agat_parser` | no | derived from `ensembl_genes_repo` | Path to `parse_agat.py` |
| `--run_agat_metrics` | no | `true` | Enable the AGAT metrics branch |
| `--outdir` | no | `./results` | Output directory |

## Running

From the repository root:

```bash
nextflow run pipelines/qc/main.nf \
  --gff_csv /path/to/gff_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --outdir results  \
  -process.executor slurm

```

If the parser or feature-levels YAML are not present under the supplied `ensembl-genes` checkout, pass them explicitly:

```bash
nextflow run pipelines/qc/main.nf \
  --gff_csv /path/to/gff_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --agat_parser /path/to/parse_agat.py \
  --feature_levels /path/to/feature_levels.yaml \
  --outdir results
```

## Outputs

Published outputs are written to:

```text
<outdir>/qc/agat/
```

        val  feature_levels_yaml
Per sample, the pipeline currently produces:

- `<sample>_agat_stats.txt`: raw AGAT statistics output
- `<sample>_agat_stats_genebuild.csv`: parsed metrics in CSV format

## Validation

The workflow stops early if any of the following are missing or invalid:

- `--gff_csv`
- `--ensembl_genes_repo`
- derived or explicit `feature_levels.yaml`
- derived or explicit `parse_agat.py`

## Implementation Notes

- The workflow entrypoint is `pipelines/qc/main.nf`.
- AGAT statistics are implemented in `modules/agat/run_agat_stats.nf`.
- AGAT report parsing is implemented in `modules/agat/parse_agat.nf`.
- The AGAT subworkflow is defined in `subworkflows/agat/agat_stats.nf`.
- `AGAT_RUN_STATS` and `AGAT_PARSE` are configured with `maxForks = 8` in `pipelines/qc/nextflow.config`.

## Current Scope

This README reflects the current implementation in this branch. The QC pipeline currently exposes the AGAT metrics path only; additional QC modules can be documented here as they are added.
