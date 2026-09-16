---
name: ensembl-nf-implement-pipeline
description: Implement or extend a pipeline in the ensembl-genes-nf repository. Use when adding a pipeline, workflow entry point, pipeline-local configuration, schema parameter, or cross-component feature while preserving this repository's DSL2 architecture and testing conventions.
---

# Implement an Ensembl Nextflow Pipeline

Inspect the target pipeline and its nearest working analogue before editing. Treat `pipelines/example/` as the learning template and an existing pipeline as the implementation precedent; do not copy a template's unfinished TODOs into production code.

## Preserve the architecture

- Keep pipeline-owned code beneath `pipelines/<name>/`: `main.nf`, `nextflow.config`, modules, subworkflows, resources, and tests/examples as needed.
- Use DSL2 and let `main.nf` compose named subworkflows. Keep tool processes out of the entry point except for truly pipeline-specific minimal work.
- Prefer one primary bioinformatics tool per module when that improves reuse and clarity. Put orchestration and channel transformations in subworkflows where useful.
- Pass sample data as tuples beginning with `meta`; preserve `meta` through outputs. State tuple shapes in comments at public subworkflow boundaries.
- Keep reusable modules configurable through `task.ext` and pipeline configuration. Put a pipeline-only parameter, resource override, or profile in that pipeline's `nextflow.config`; do not modify root configuration for a one-pipeline need.
- Add or update user-facing parameters in the pipeline schema/configuration together when those artifacts are used. Validate required inputs early when practical, before work is launched.

## Design gates before editing

- Read the target pipeline's recent history with `git log` and identify every
  execution mode, optional resource, and phase boundary before changing code.
  In the long-read TAMA history, `c4a9ed3` and `82127af` added modes to
  `main.nf`, while `2a6fd19` moved processes behind subworkflows without
  removing the entry-point state machine. Treat that sequence as a warning:
  shortening `main.nf` is not sufficient if phase boundaries and mode contracts
  remain implicit.
- Keep `main.nf` as a composition boundary: process modules belong behind named
  subworkflows. A terminal software-version aggregation process may remain a
  direct entry-point call when that matches an existing repository precedent.
- Treat `nextflow_schema.json` as an executable contract. The entry workflow
  must call `validateParameters()`, use the nf-schema-supported JSON Schema
  dialect, describe inherited reporting parameters, and express types, enums,
  file existence, and conditional requirements where possible. Keep only
  genuinely cross-field/runtime rules in a small custom validator.
- For each conditional process or optional database/index, design an explicit
  mode matrix. Initialize absent outputs and versions with `channel.empty()`;
  never reference an output from a process that may not be invoked. Make
  singleton files value/broadcast channels before pairing them with sample
  streams.
- Separate inventory/approval from production processing with named workflows
  and explicit emitted contracts. Do not make `return` the only indication that
  a phase stopped.

## Work deliberately

1. Read the affected `main.nf`, `nextflow.config`, direct modules/subworkflows, and the relevant sections of `docs/NEXTFLOW_REQUIREMENTS.md`, `docs/template/MODULES.md`, `docs/template/PATTERNS.md`, and `docs/template/CONFIGURATION.md`.
2. Trace changed input and output channels end-to-end. Name emitted channels for their content, not their position.
3. Reuse resource labels from `config/resources.config`, explicit tool versions, and an appropriate container/Conda definition where applicable.
4. Keep outputs under `${params.outdir}` and preserve standard execution reporting through the inherited root configuration.
5. Add a small, representative stub-mode path for changed orchestration when useful. Do not add generated results, `.nextflow*`, `work/`, or container caches to version control.

## Finish with evidence

Run the narrowest relevant commands first:

```bash
nextflow lint -o concise pipelines/<pipeline>
nextflow config pipelines/<pipeline>/main.nf
nextflow run pipelines/<pipeline>/main.nf -stub-run -profile test --outdir /tmp/<pipeline>-stub
```

Use the pipeline's established test command when it exists. If an external tool, image, reference, or test input prevents execution, report the exact command, blocker, and the unverified surface; never claim a real run passed from a parse-only check.

For every new or changed optional mode, run a stub/configuration matrix that
covers the option supplied and omitted, plus the inventory-only path when one
exists. Include at least one intentionally invalid parameter combination and
verify that nf-schema rejects it before any process is submitted.

## Production safety gates

For maintained production pipelines:

- Require the documented minimum Nextflow version; this repository currently requires 26.04.6 or newer.
- Enable timeline, report, and trace outputs for production runs; DAG output is optional.
- Use `set -euo pipefail` or an equivalent global `process.shell` for every real process script.
- Never ignore failures in required production paths. Retry only acquisition or infrastructure-sensitive work; terminate on validation, alignment, TAMA, and model-QC failures.
- Put diagnostic messages before assertions for required inputs and outputs.
- Make helper-script invocation reproducible by either resolving an executable from a configured pipeline path or intentionally staging it as a declared input.
- A stub run proves wiring only; completion also requires focused failure tests and any requested real container/HPC execution.
