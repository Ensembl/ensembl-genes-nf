---
name: ensembl-nf-implement-pipeline
description: Implement or extend a pipeline in the ensembl-genes-nf repository. Use when adding a pipeline, workflow entry point, pipeline-local configuration, schema parameter, or cross-component feature while preserving this repository's DSL2 architecture and testing conventions.
---

# Implement an Ensembl Nextflow Pipeline

Inspect the target pipeline and its nearest working analogue before editing. Treat `pipelines/example/` as the learning template and an existing pipeline as the implementation precedent; do not copy a template's unfinished TODOs into production code.

## Preserve the architecture

- Keep pipeline-owned code beneath `pipelines/<name>/`: `main.nf`, `nextflow.config`, modules, subworkflows, resources, and tests/examples as needed.
- Use DSL2 and let `main.nf` compose named subworkflows. Keep tool processes out of the entry point except for truly pipeline-specific minimal work.
- Use one primary bioinformatics tool per module. Put orchestration and channel transformations in subworkflows.
- Pass sample data as tuples beginning with `meta`; preserve `meta` through outputs. State tuple shapes in comments at public subworkflow boundaries.
- Keep reusable modules configurable through `task.ext` and pipeline configuration. Put a pipeline-only parameter, resource override, or profile in that pipeline's `nextflow.config`; do not modify root configuration for a one-pipeline need.
- Add or update user-facing parameters in the pipeline schema/configuration together. Validate required inputs early, before work is launched.

## Work deliberately

1. Read the affected `main.nf`, `nextflow.config`, direct modules/subworkflows, and the relevant sections of `docs/NEXTFLOW_REQUIREMENTS.md`, `docs/template/MODULES.md`, `docs/template/PATTERNS.md`, and `docs/template/CONFIGURATION.md`.
2. Trace every changed input and output channel end-to-end. Name emitted channels for their content, not their position.
3. Reuse resource labels from `config/resources.config`, explicit tool versions, and an appropriate container/conda definition.
4. Keep outputs under `${params.outdir}` and preserve standard execution reporting through the inherited root configuration.
5. Add a small, representative stub-mode path for changed orchestration. Do not add generated results, `.nextflow*`, `work/`, or container caches to version control.

## Finish with evidence

Run the narrowest relevant commands first:

```bash
nextflow lint -o concise pipelines/<pipeline>
nextflow config pipelines/<pipeline>/main.nf
nextflow run pipelines/<pipeline>/main.nf -stub-run -profile test --outdir /tmp/<pipeline>-stub
```

Use the pipeline's established test command when it exists. If an external tool, image, reference, or test input prevents execution, report the exact command, blocker, and the unverified surface; never claim a real run passed from a parse-only check.
