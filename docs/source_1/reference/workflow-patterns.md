# Workflow Patterns

The repository uses a few recurring Nextflow DSL2 patterns. Keeping these consistent makes it easier to move between pipelines.

## Public Entry Point

Use `main.nf` as the public entry point for each pipeline. Keep parsing, validation, and top-level orchestration there; place reusable logic in modules and subworkflows.

## Metadata Maps

Pass sample context as a `meta` map in tuple channels:

```groovy
tuple val(meta), path(input_file)
```

Common keys include `id`, `gca`, `species_name`, `study_accession`, and tool-specific fields. Preserve the map through modules so output files can be grouped and named consistently.

## Parameter Validation

Production pipelines should include `nf-schema` validation:

```groovy
include { validateParameters } from 'plugin/nf-schema'
validateParameters()
```

Use `nextflow_schema.json` as the source of truth for parameter names, defaults, required values, and descriptions.

## Channel Composition

Prefer small modules with explicit `emit:` names, then compose them in subworkflows:

```groovy
emit:
genome_bam = SAMTOOLS_INDEX_GENOME.out.bam_and_bai
logs = STAR_ALIGN.out.log
```

Named emits make downstream wiring and documentation much easier to understand.

## Publishing

Use `publishDir` for user-facing outputs and keep intermediate files inside the work directory unless they are needed for inspection, reporting, or reuse. Published paths should be predictable and usually live under `params.outdir`.

## Version Tracking

Modules should emit a `versions.yml` file where practical. This gives downstream reports and release audits a consistent place to find tool versions.
