# Nextflow requirements

This repository targets Nextflow 25.10.2 or newer, Java 17 or newer, and
strict/v2 syntax. Nextflow 26.04 enables the v2 parser by default; on
Nextflow 25.10.2 through 25.x set `NXF_SYNTAX_PARSER=v2`.

The following is the recommended baseline for maintained pipelines. It is a
starting point rather than a checklist that every pipeline must follow without
exception. If a pipeline needs a different choice, document the reason in its
README and keep the test and user-facing behavior clear.

The usual defaults are:

1. Use DSL2 with an explicit entry `workflow`.
2. Put executable statements inside a process, workflow, or function.
3. Validate parameters inside the entry workflow, where practical.
4. Prefer `channel` over the deprecated `Channel` factory.
5. Prefer `script:` for process commands and declare process inputs explicitly.
6. Keep high-level branching in workflows and tool execution in modules when
   that makes the code easier to reuse.
7. Pass reference files and helper scripts as process inputs when they are part
   of the task's data or need to be staged reproducibly.
8. Provide a stub implementation for workflows that need quick structural
   testing, with outputs that match the real process.
9. Keep the schema, README, defaults, and test profile broadly consistent.

## Checks

For the template example:

```bash
cd pipelines/example
nextflow lint -o concise .
nextflow run main.nf -stub-run -profile test
```

For Nextflow 25.10.2 through 25.x:

```bash
cd pipelines/example
NXF_SYNTAX_PARSER=v2 nextflow lint -o concise .
NXF_SYNTAX_PARSER=v2 nextflow run main.nf -stub-run -profile test
```

The repository also contains established production pipelines. Run a
repository-wide lint as a separate migration check when changing those
pipelines; existing production code may not yet satisfy strict syntax.

The root configuration selects Singularity and disables Docker and Conda. The
Ensembl HPC cache location is not part of the portable default; apply
`config/singularity.config` with `-c` on the cluster.

```bash
cd ../..
nextflow lint -o concise .
```

For a new pipeline, a small test profile and deterministic input under the
pipeline directory are useful defaults. Avoid production data and
host-specific paths in tests where possible.
