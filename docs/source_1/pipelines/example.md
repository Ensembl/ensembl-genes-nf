# Example Pipeline

The example pipeline is a small runnable workflow for learning the repository patterns. It demonstrates single-process workflows, subworkflow composition, module outputs, and stub-mode development.

## Run It

```bash
cd pipelines/example
nextflow run main.nf -stub --input example.txt --outdir results
```

For the simplest possible workflow, run:

```bash
nextflow run workflows/simple_workflow.nf -stub --outdir results
```

## What It Demonstrates

| Area | What to inspect |
| --- | --- |
| Simple workflow | `pipelines/example/workflows/simple_workflow.nf` |
| Subworkflow composition | `pipelines/example/main.nf` and `pipelines/example/subworkflows/` |
| Process modules | `pipelines/example/modules/` |
| Schema validation | `pipelines/example/nextflow_schema.json` |

Use this pipeline as the first place to test new conventions before applying them to a production workflow.

## Parameters

```{toctree}
:maxdepth: 1

../generated/example-parameters
```
