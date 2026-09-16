---
name: ensembl-nf-build-component
description: Create or modify a conformant Nextflow DSL2 module or subworkflow in ensembl-genes-nf. Use for process wrappers, channel contracts, module reuse, subworkflow composition, containers, resource labels, versions files, and stubs.
---

# Build a Conformant Nextflow Component

Read `docs/template/MODULES.md`, `modules/reference_example.nf`, and the closest real component before changing code. Prefer the real pipeline's established interface where it is more specific than the template.

## Recommended module defaults

- Give a module one primary tool responsibility when practical. Split unrelated tools into separate modules when that improves reuse or clarity.
- Accept and emit `tuple val(meta), path(...)` for sample-associated files; retain the same `meta` map without mutation. Create a new map with `meta + [...]` only when metadata must change.
- Set a meaningful `tag`, resource `label`, and version-pinned container or Conda environment when the process needs them. Add `publishDir` when outputs are user-facing and the module owns publication policy.
- Emit named result channels. Add `path "versions.yml", emit: versions` when version reporting is part of the pipeline's reproducibility approach.
- Define optional arguments and output prefixes in `task.ext` when configuration needs to vary between callers. Use `task.ext.when` when conditional execution is best controlled externally.
- Make the `stub` create every declared output when stub testing is provided, using the same naming rules as the real script.

## Subworkflow defaults

- Declare `take:`, `main:`, and `emit:` explicitly for reusable named subworkflows. Keep public input/output tuple shapes in comments.
- Use descriptive channel names. Apply `join`, `groupTuple`, `collect`, `branch`, or `ifEmpty` only after checking the expected cardinality and metadata keys.
- Expose useful result and version channels; avoid hiding outputs required by later workflow stages.
- Keep process resource settings and tool arguments in configuration, not hardcoded in the subworkflow.

## Verification

Parse the owning pipeline and exercise the smallest stub-mode workflow that reaches the component when a stub path exists. Inspect output tuple names and generated version reports when they are part of the component. Keep fixes scoped: do not reformat unrelated modules or replace a working local convention with an abstract template.

## Failure and reproducibility gates

For maintained production components, also verify:

- The process has one primary tool or one inseparable helper responsibility; split independent acquisition, validation, conversion, alignment, sorting, and QC concerns.
- Required inputs and declared outputs have explicit user-facing failure messages before checks.
- Real scripts run with `-euo pipefail`, and a pipe is tested for an upstream failure that a downstream tool could mask.
- `errorStrategy 'ignore'` is absent from required stages; retries are limited to named acquisition/infrastructure processes.
- Stub output names and tuple arity exactly match the real contract.
- Reports, checksums, audits, and versions are declared outputs, not incidental work-directory files.
- Structured data is parsed structurally rather than with fragile regular-expression `sed` extraction.
- At least one generated `.command.sh` is inspected for quoting, escaping, shell options, and helper-script staging.
