---
name: ensembl-nf-verify-change
description: Review and validate a change in ensembl-genes-nf before handoff. Use for Nextflow syntax/configuration checks, stub runs, channel-contract review, schema/configuration consistency, output hygiene, and focused repository change review.
---

# Verify an Ensembl Nextflow Change

Establish the changed surface with `git diff --check`, `git status --short`, and a targeted diff. Treat unrelated dirty files as user work: do not stage, delete, or edit them.

Run these deterministic structural checks before deeper review:

```bash
test ! -f pipelines/<pipeline>/nextflow_schema.json || rg -q 'validateParameters\(\)' pipelines/<pipeline>/main.nf
test ! -f pipelines/<pipeline>/nextflow_schema.json || rg -q 'https://json-schema.org/draft/2020-12/schema' pipelines/<pipeline>/nextflow_schema.json
test -z "$(rg -n 'include.*from .*(\.\./)?modules/' pipelines/<pipeline>/main.nf | rg -v 'COLLECT' || true)"
```

The final check intentionally permits no direct process-module imports in the
entry point. If a repository-standard terminal aggregation process is a
documented exception, state that exception explicitly in the review.

## Review against repository rules

- Check that DSL2 entry points and subworkflows are structured clearly, and that module boundaries are sensible for the tools they wrap.
- Trace each changed tuple from producer to consumer: `meta` is retained, file arity matches, and `emit` names match their consumers.
- Check applicable process conventions: resource labels, version-pinned container/Conda definitions, `task.ext` defaults, version reporting, and faithful stub blocks.
- Confirm pipeline-specific params, tool arguments, resource overrides, and profiles live in `pipelines/<name>/nextflow.config`; root configuration remains shared only. Keep HPC-specific cache paths in an explicit config such as `config/singularity.config`.
- Confirm declared output filenames, `publishDir` patterns, and stub outputs agree where those features are used. Generated results, `work/`, `.nextflow*`, caches, and local reference paths must stay out of the diff.
- Inspect recent commits when reviewing a refactor. Identify newly added modes
  and verify that each has a named workflow boundary, schema coverage, and a
  corresponding stub or failure test.

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

For every conditional process/resource, run both sides of its mode matrix. At
minimum this means: optional resource supplied and omitted; optional branch
enabled and disabled; inventory-only and production modes when applicable; and
one invalid combination expected to fail through nf-schema. Confirm that the
disabled branch produces no task and no undefined output-channel access.

## Required production-path audit

Before handoff, explicitly audit the changed pipeline for:

- Nextflow minimum version >= 26.04.6 and enabled timeline/report/trace outputs.
- No `errorStrategy 'ignore'` in required production processes; retries restricted to acquisition or infrastructure-sensitive tasks.
- `-euo pipefail` in generated real scripts, with a generated `.command.sh` inspected for shell correctness.
- Explicit failure messages for corrupt gzip, malformed FASTQ, checksum mismatch, missing sidecars, empty approved manifests, empty TAMA inputs/outputs, and invalid model rows when applicable.
- Public channel shapes traced producer-to-consumer, including preserved `meta` fields and file arity.
- Conditional output references are branch-safe, singleton inputs are reusable
  across all intended consumers, and every emitted version channel is included
  in the final software-version report.
- Focused failure tests for every changed validation/tool boundary, plus a real container/Slurm run when acceptance criteria require it.

Do not report a pipeline as complete if only lint, config parsing, or stub execution has passed.
