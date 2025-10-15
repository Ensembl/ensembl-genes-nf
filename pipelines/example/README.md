# Example Pipeline

Minimal Nextflow pipeline demonstrating subworkflow patterns.

## Quick Start

```bash
# Run the main workflow (2 tools + combine + count)
nextflow run main.nf -stub --outdir results
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
├── advanced_entrypoints/         # Advanced: Entry point system
└── docs/                         # All documentation
```

## Examples

### Simple Workflow (workflows/simple_workflow.nf)
A single process workflow - demonstrates the basics.

```bash
nextflow run workflows/simple_workflow.nf -stub --outdir results
```

**Output**: `tool_a/` (1 tool, 2 samples)

### Main Workflow (main.nf)
Two subworkflows chained together (also available as `workflows/subworkflow_example.nf`).

```bash
nextflow run main.nf -stub --outdir results
```

**Output**: `tool_a/`, `tool_b/`, `combined/`, `line_counts/` (2 tools + processing)

### Advanced Entry Points (advanced_entrypoints/)
Dynamic entry point system with automatic workflow resumption.

```bash
cd advanced_entrypoints
nextflow run main.nf --help
nextflow run main.nf -stub --outdir results
```

**Output**: All 4 tools + combine + count with automatic entry point detection

## Where to Start

**New to Nextflow?** Start with `workflows/simple_workflow.nf` to understand the basics.

**Familiar with workflows?** Explore `main.nf` to see subworkflow patterns.

**Building complex pipelines?** Check out `advanced_entrypoints/` for the entry point system.

## Documentation

All documentation is in **[docs/](docs/)**:
- [docs/QUICK_START.md](docs/QUICK_START.md) - Complete guide
- [docs/PATTERNS.md](docs/PATTERNS.md) - Common patterns
- [docs/INDEX.md](docs/INDEX.md) - Full index

## Parameters

- `--outdir` - Output directory (default: `'results'`)
