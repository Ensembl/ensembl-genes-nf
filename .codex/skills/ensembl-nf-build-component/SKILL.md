---
name: ensembl-nf-build-component
description: Create or modify a conformant Nextflow DSL2 module or subworkflow in ensembl-genes-nf. Use for process wrappers, channel contracts, module reuse, subworkflow composition, containers, resource labels, versions files, and stubs.
---

# Build a Conformant Nextflow Component

Read `docs/template/MODULES.md`, `modules/reference_example.nf`, and the closest real component before changing code. Prefer the real pipeline's established interface where it is more specific than the template.

## Module contract

- Give a module one primary tool responsibility. Split unrelated tools into separate modules and compose them in a subworkflow.
- Accept and emit `tuple val(meta), path(...)` for sample-associated files; retain the same `meta` map without mutation. Create a new map with `meta + [...]` only when metadata must change.
- Set a meaningful `tag`, resource `label`, version-pinned container or conda environment, and an explicit `publishDir` under `${params.outdir}` when outputs are user-facing.
- Emit named result channels and `path "versions.yml", emit: versions`. Capture the real tool version in `script` and a fixed representative value in `stub`.
- Define optional arguments and output prefix in `task.ext` with safe defaults. Use `task.ext.when` for externally configured conditional execution instead of embedding pipeline policy in the module.
- Make the `stub` create every declared output with the same naming rules as the real script.

## Subworkflow contract

- Declare `take:`, `main:`, and `emit:` explicitly. Keep public input/output tuple shapes in comments.
- Use descriptive channel names. Apply `join`, `groupTuple`, `collect`, `branch`, or `ifEmpty` only after checking the expected cardinality and metadata keys.
- Expose useful result and version channels; avoid hiding outputs required by later workflow stages.
- Keep process resource settings and tool arguments in configuration, not hardcoded in the subworkflow.

## Verification

Parse the owning pipeline and exercise the smallest stub-mode workflow that reaches the component. Inspect the output tuple names and generated `versions.yml` before declaring the work complete. Keep fixes scoped: do not reformat unrelated modules or replace a working local convention with an abstract template.
