# Getting Started

## Requirements

The pipelines are written in Nextflow DSL2 and are designed to run with containerised tools. A typical production run needs:

- Nextflow.
- Singularity or another supported container engine, depending on the profile.
- Access to the input files, reference data, and any Ensembl/EBI file systems used by the selected pipeline.
- A pipeline-specific `nextflow_schema.json` file kept in sync with the accepted parameters.

## Choose a Pipeline

Start from the pipeline directory rather than the repository root. Each pipeline has its own `main.nf`, `nextflow.config`, schema, modules, and sometimes extra assets.

```bash
nextflow run pipelines/repeat/main.nf \
  -c pipelines/repeat/nextflow.config \
  --csvFile samples.csv \
  --outdir results/repeat \
  --generate_lib \
  --run_repeatmasker
```

The exact parameters differ by pipeline. The generated parameter pages in this documentation are built from the schema files so that the docs and validation rules stay close together.

## Validate Inputs Early

Most production pipelines use `nf-schema` validation. When you add or rename a parameter, update the pipeline schema at the same time as the code. For CSV inputs, keep the separate CSV schema near the pipeline assets when available.

## Outputs

The shared root `nextflow.config` enables timeline, report, trace, and DAG outputs under `params.tracedir`. Pipeline-specific outputs are usually published under `params.outdir` with per-sample or per-tool subdirectories.

## Local Documentation Build

Install the documentation dependencies, regenerate parameter pages, and build HTML:

```bash
python -m pip install -r docs/requirements.txt
python docs/scripts/render_schema_docs.py
sphinx-build -b html docs/source docs/build/html
```

Open `docs/build/html/index.html` in a browser after the build completes.
