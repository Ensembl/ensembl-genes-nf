# Configuration

## Base Configuration

The root `nextflow.config` defines shared defaults for all pipelines:

- Repository manifest metadata.
- Default `params.outdir` and `params.tracedir`.
- Timeline, report, trace, and DAG generation.
- Singularity defaults.
- Shared resource configuration from `config/resources.config` and `config/shared-module.config`.

These settings provide a consistent starting point. Pipeline-specific configuration files can override them when the workflow needs different execution profiles, containers, parameters, or resource labels.

## Pipeline Parameters

Each pipeline should keep a `nextflow_schema.json` file beside its `main.nf`. The schema is used by `nf-schema` for runtime validation and by this documentation site to generate the parameter reference.

When changing a parameter:

1. Update the default in `nextflow.config` when there is one.
2. Update `nextflow_schema.json` with type, description, default, and required status.
3. Regenerate docs with `python docs/scripts/render_schema_docs.py`.
4. Add or update an example command in the relevant pipeline page when the user-facing behaviour changes.

## Profiles

Use profiles for execution environments rather than changing workflow logic. Existing pipelines define combinations such as `slurm`, `local`, `conda`, `docker`, and `singularity`.

```bash
nextflow run pipelines/riboseq/main.nf \
  -c pipelines/riboseq/nextflow.config \
  -profile slurm \
  --sample_sheet samples.csv \
  --star_index /path/to/star_index \
  --gtf annotation.gtf \
  --fasta genome.fa \
  --chrom_sizes_file genome.chrom.sizes
```

## Containers and Caches

The repository defaults to Singularity for HPC execution. Some modules use BioContainers, GHCR images, or Wave/ORAS references. Keep cache paths explicit in configuration so large images do not end up in user home directories by accident.

## Reporting

The shared reports are written under `params.tracedir`:

| Report | Default file |
| --- | --- |
| Timeline | `results/pipeline_info/execution_timeline.html` |
| Execution report | `results/pipeline_info/execution_report.html` |
| Trace | `results/pipeline_info/execution_trace.txt` |
| DAG | `results/pipeline_info/pipeline_dag.svg` |
