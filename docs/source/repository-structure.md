# Repository Structure

The repository is arranged as a collection of pipelines plus shared building blocks.

```text
ensembl-genes-nf/
├── config/                 Shared process labels and module configuration
├── docs/                   Sphinx documentation source and older Markdown guides
├── modules/                Shared module examples and cross-pipeline helpers
├── pipelines/
│   ├── example/            Minimal pattern-focused workflow
│   ├── repeat/             Repeat annotation workflow
│   ├── riboseq/            Ribosome profiling workflow
│   └── translon-consensus/ Consensus reporting workflow
├── subworkflows/           Shared subworkflow examples
├── test/                   Small module or workflow test fixtures
└── nextflow.config         Base configuration inherited by pipeline configs
```

## Pipeline Layout

Most pipeline directories follow this shape:

```text
pipelines/<pipeline>/
├── main.nf
├── nextflow.config
├── nextflow_schema.json
├── modules/
├── subworkflows/ or workflows/
├── assets/
└── bin/ or scripts/
```

Use `main.nf` as the public entry point. Put reusable process definitions in `modules/`, larger DSL2 composition in `subworkflows/` or `workflows/`, and non-code templates or schemas in `assets/`.

## Shared Configuration

The repository root `nextflow.config` provides common reporting, base resources, and container defaults. Pipeline configs include it with:

```groovy
includeConfig '../../nextflow.config'
```

Pipeline configs should override only the settings that are specific to that workflow.
