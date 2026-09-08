---
name: ensembl-nf-verify-change
description: Review and validate a change in ensembl-genes-nf before handoff. Use for Nextflow syntax/configuration checks, stub runs, channel-contract review, schema/configuration consistency, output hygiene, and focused repository change review.
---

# Verify an Ensembl Nextflow Change

Establish the changed surface with `git diff --check`, `git status --short`, and a targeted diff. Treat unrelated dirty files as user work: do not stage, delete, or edit them.

## Review against repository rules

- Confirm DSL2 entry points compose subworkflows and that a module has one primary tool responsibility.
- Trace each changed tuple from producer to consumer: `meta` is retained, file arity matches, and `emit` names match their consumers.
- Confirm process settings use a resource label, version-pinned container/conda, `task.ext` defaults where applicable, a `versions.yml` output, and a faithful stub block.
- Confirm pipeline-specific params, tool arguments, resource overrides, and profiles live in `pipelines/<name>/nextflow.config`; root configuration remains shared only. Keep HPC-specific cache paths in an explicit config such as `config/singularity.config`.
- Confirm declared output filenames, `publishDir` patterns, and stub outputs agree. Generated results, `work/`, `.nextflow*`, caches, and local reference paths must stay out of the diff.

## Test in layers

Run checks from cheapest to most representative:

```bash
git diff --check
nextflow lint -o concise pipelines/<pipeline>
nextflow config pipelines/<pipeline>/main.nf
nextflow run pipelines/<pipeline>/main.nf -stub-run -profile test --outdir /tmp/<pipeline>-stub
```

Use a targeted existing test where one exists, for example:

```bash
cd pipelines/riboseq && ./test_unique_reads.sh
```

Use an isolated temporary output directory for ad-hoc runs. State exactly which layers passed and which were not attempted, including missing binaries, containers, references, or input data. A syntax/configuration check is not an execution result.
