# Nextflow conventions and validation

The example pipeline is developed and tested with Nextflow 26.04 or newer,
Java 17 or newer, and strict syntax. Nextflow 26.04 enables the v2 parser by
default. Nextflow 25.10.2 or newer may also work when `NXF_SYNTAX_PARSER=v2`
is enabled.

The following are recommended defaults for new and maintained pipelines. They
are a starting point, not a universal checklist. If a pipeline needs a
different choice, document the reason where it will help future maintainers and
keep the test and user-facing behavior clear.

The usual defaults are:

1. Use DSL2 with an explicit entry `workflow` for new pipelines.
2. Put executable statements inside a process, workflow, or function where
   possible.
3. Validate parameters inside the entry workflow when that makes ownership
   clear.
4. Prefer `channel` over the older `Channel` factory in new code.
5. Prefer `script:` for process commands and declare process inputs explicitly
   when applicable.
6. Keep high-level branching in workflows and tool execution in modules when
   that makes the code easier to reuse.
7. Pass reference files and helper scripts as process inputs when they are part
   of the task's data or need to be staged reproducibly.
8. Provide a stub implementation when quick structural testing is useful, with
   outputs that match the real process.
9. Keep the schema, README, defaults, and test profile broadly consistent when
   those artifacts are present.

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
