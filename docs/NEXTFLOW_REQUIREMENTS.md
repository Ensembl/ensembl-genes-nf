# Nextflow requirements

This repository targets Nextflow 26.04 or newer, Java 17 or newer, and strict
syntax. Nextflow 26.04 enables the strict parser by default.

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

From the repository root:

```bash
nextflow lint -o concise .
```

From the example pipeline:

```bash
cd pipelines/example
nextflow lint -o concise .
nextflow run main.nf -stub-run -profile test
```

For a new pipeline, a small test profile and deterministic input under the
pipeline directory are useful defaults. Avoid production data and
host-specific paths in tests where possible.
