# Example Pipeline

Minimal Nextflow pipeline demonstrating subworkflow patterns.

## Quick Start

```bash
# Run the main workflow (2 tools + combine + count)
nextflow run main.nf -stub-run -profile test
```

This runs two tools in parallel, combines their outputs, and counts lines.

## Structure

```
.
├── main.nf                       # Subworkflow example (2 tools)
├── workflows/
│   ├── simple_workflow.nf       # Simplest: 1 process only
│   └── subworkflow_example.nf   # Same as main.nf (for reference)
├── subworkflows/                 # 2 reusable subworkflow examples
├── modules/                      # 7 example processes
├── assets/                       # Deterministic test inputs
└── docs/                         # All documentation
```

## Examples

### Simple Workflow (workflows/simple_workflow.nf)
A single process workflow - demonstrates the basics.

```bash
nextflow run workflows/simple_workflow.nf -stub-run --outdir results
```

**Output**: `tool_a/` (1 tool, 2 samples)

### Main Workflow (main.nf)
Two subworkflows chained together (also available as `workflows/subworkflow_example.nf`).

```bash
nextflow run main.nf -stub-run -profile test
```

**Output**: `tool_a/`, `tool_b/`, `combined/`, `line_counts/` (2 tools + processing)

## Where to Start

**New to Nextflow?** Start with `workflows/simple_workflow.nf` to understand the basics.

**Familiar with workflows?** Explore `main.nf` to see subworkflow patterns.

**Building complex pipelines?** Use the workflow/subworkflow/module separation
shown in `main.nf` and the `workflows/`, `subworkflows/`, and `modules/`
directories.

## Documentation

All documentation is in **[docs/](docs/)**:
- [docs/QUICK_START.md](docs/QUICK_START.md) - Complete guide
- [docs/PATTERNS.md](docs/PATTERNS.md) - Common patterns
- [docs/INDEX.md](docs/INDEX.md) - Full index

## Parameters

- `--input` - Input file (default: `assets/input.txt`)
- `--outdir` - Output directory (default: `./results`)

The schema, defaults, README, and workflow are intentionally kept in sync so
this directory can be copied as a starting point for a real pipeline.

## Why `main.nf`, `workflows/`, and `subworkflows/` are separate

The separation keeps command-line parameters at the boundary of the pipeline.
`main.nf` is responsible for validating named parameters such as
`params.input`, resolving files and defaults, and converting them into the
channels and values that the pipeline actually needs. The workflow layer then
decides which stages run and how those inputs are connected.

For reusable subworkflows, it is usually clearer not to reach directly into
named parameters. A subworkflow can instead receive explicit inputs in its
`take:` block, such as `samples_ch`,
`reference`, or `mode`, and expose explicit outputs in its `emit:` block. This
means a subworkflow does not silently depend on a parameter name, a particular
schema, or a particular entrypoint. A different pipeline can reuse the same
subworkflow with a different parameter schema, or call it with a channel built
from a manifest, database, or another upstream stage.

In practice, the dependency flow is:

```text
named CLI/config parameters
        │
        ▼
main.nf: validate, resolve, normalize
        │  explicit channels and values
        ▼
workflows/: choose stages and branches
        │  explicit subworkflow inputs
        ▼
subworkflows/: compose reusable analysis units
        │  explicit module inputs
        ▼
modules/: execute one tool or transformation
```

For example, a QC entrypoint may turn `--mode combined` and an annotation
manifest into `annotation_ch`, `database`, and `data_file_path`. The QC
subworkflow can then take those values explicitly without knowing whether they
came from command-line parameters, a profile, or another workflow. This keeps
parameter dependencies manageable at the edge instead of imposing one
pipeline's parameter names on every reusable component. Direct parameter access
can still be appropriate for genuinely global settings; the useful question is
whether the dependency is intentional and easy to see.
