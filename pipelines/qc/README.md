# QC Pipeline

This pipeline runs annotation QC on one or more samples using AGAT statistics and, optionally, InterProScan-derived metrics.

## What It Does

For each row in the input CSV, the pipeline:

1. runs `agat_sp_statistics.pl` on the GFF3 file
2. parses the AGAT output into a `*_agat_stats_genebuild.csv` file
3. optionally runs InterProScan on the protein FASTA
4. optionally parses the InterProScan output into an InterPro metrics TSV

## Requirements

- Nextflow 25.x
- Singularity available on the execution host
- InterProScan 5.78-109.0 data available on the execution host when the InterProScan branch is enabled
- A local checkout of `ensembl-genes` containing:
    - `src/python/ensembl/genes/annotation_qc/parsers/parse_agat.py`
    - `src/python/ensembl/genes/annotation_qc/config/feature_levels.yaml`
    - `src/python/ensembl/genes/annotation_qc/parsers/interpro.py` when `--run_interproscan` is enabled

## Inputs

### Sample Sheet

The required input is `--input_csv`, a CSV file with a header. Include the columns needed by the enabled branches:

```csv
sample,gff3,protein
sample_1,/path/to/sample_1.gff3,/path/to/sample_1.faa
sample_2,/path/to/sample_2.gff3,/path/to/sample_2.faa
```

- `sample`: sample identifier used for logging and task tags
- `gff3`: path to the GFF3 file used by AGAT
- `protein`: path to the protein FASTA used by InterProScan

If `--run_agat_metrics` is enabled, `gff3` must be present for every row.
If `--run_interproscan` is enabled, `protein` must be present for every row.

## Parameters

| Parameter              | Required | Default                                                                   | Description                               |
| ---------------------- | -------- | ------------------------------------------------------------------------- | ----------------------------------------- |
| `--input_csv`          | yes      | none                                                                      | CSV describing the samples and file paths |
| `--ensembl_genes_repo` | yes      | none                                                                      | Path to a local `ensembl-genes` checkout  |
| `--feature_levels`     | no       | derived from `ensembl_genes_repo`                                         | Path to `feature_levels.yaml`             |
| `--agat_parser`        | no       | derived from `ensembl_genes_repo`                                         | Path to `parse_agat.py`                   |
| `--run_agat_metrics`   | no       | `false`                                                                   | Enable the AGAT metrics branch            |
| `--run_interproscan`   | no       | `false`                                                                   | Enable the InterProScan branch            |
| `--database`           | no       | `Pfam`                                                                    | InterProScan application database         |
| `--data_file_path`     | no       | none                                                                      | Path to the InterProScan 5 data directory |
| `--interpro_parser`    | no       | derived from `ensembl_genes_repo`                                         | Path to `interpro.py`                     |
| `--outdir`             | no       | `./results`                                                               | Output directory                          |

## Running

```bash
nextflow run pipelines/qc/main.nf \
  --input_csv /path/to/qc_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --outdir results \
  -profile slurm
```

To enable AGAT:

```bash
nextflow run pipelines/qc/main.nf \
  --input_csv /path/to/qc_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --run_agat_metrics \
  --outdir results \
  -profile slurm
```

To enable InterProScan:

```bash
nextflow run pipelines/qc/main.nf \
  --input_csv /path/to/qc_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --run_interproscan \
  --database Pfam \
  --data_file_path /path/to/interproscan/data \
  --outdir results \
  -profile slurm
```

The `slurm` profile enables Singularity and runs the QC processes through Slurm. The InterProScan container is pulled from the Docker image `interpro/interproscan:5.78-109.0`; Docker itself is not required on the cluster.

If the parsers or feature-level YAML are not present under the supplied `ensembl-genes` checkout, pass them explicitly:

```bash
nextflow run pipelines/qc/main.nf \
  --input_csv /path/to/qc_samples.csv \
  --ensembl_genes_repo /path/to/ensembl-genes \
  --agat_parser /path/to/parse_agat.py \
  --interpro_parser /path/to/interpro.py \
  --feature_levels /path/to/feature_levels.yaml \
  --outdir results
```

## Outputs

Published outputs are written to:

```text
<outdir>/qc/agat/
<outdir>/qc/interpro/
```

Per sample, the pipeline currently produces:

- `<sample>_agat_stats.txt`: raw AGAT statistics output
- `<sample>_agat_stats_genebuild.csv`: parsed AGAT metrics
- `<sample>.tsv`: raw InterProScan TSV output
- `<sample>_hits.tsv`: parsed InterPro metrics


## Validation

`nf-schema` validates parameter types, paths, file existence, the InterProScan database allow-list, branch requirements, and the sample-sheet structure. The workflow adds only checks that require deriving paths from the supplied repository:

- derived or explicit `feature_levels.yaml`
- derived or explicit `parse_agat.py`
- derived or explicit `interpro.py` when `--run_interproscan` is enabled

## Separation of Concerns

Each layer has one job:

- The schema validates inputs: parameter types, paths, allowed database names, and sample-sheet columns.
- `main.nf` connects the input sheet to the selected branches and handles rules such as which columns are required for the enabled metrics.
- The subworkflows describe reusable analysis chains, such as running InterProScan and then parsing its output.
- The modules run one tool or parser each.

This means another pipeline can reuse the InterProScan workflow with its own input channels and output handling, without depending on this QC pipeline's command-line parameters or sample-sheet format.

## Implementation Notes

- The workflow entrypoint is `pipelines/qc/main.nf`.
- AGAT statistics are implemented in `modules/agat/run_agat_stats.nf`.
- AGAT report parsing is implemented in `modules/agat/parse_agat.nf`.
- The AGAT subworkflow is defined in `subworkflows/agat/agat_stats.nf`.
- InterProScan execution is implemented in `modules/interpro/run_interpro.nf`.
- InterPro metrics parsing is implemented in `modules/interpro/interpro_stats.nf`.
- The InterPro subworkflow is defined in `subworkflows/interpro/interproscan.nf`.

## Current Scope

This pipeline exposes AGAT metrics behind `--run_agat_metrics` and InterProScan metrics behind `--run_interproscan`.
