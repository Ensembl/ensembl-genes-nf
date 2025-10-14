# Advanced Entry Point System

This directory contains the **dynamic entry point system** - an advanced pattern for managing multi-stage workflow execution with automatic dependency resolution.

## What's Here

- **main.nf** - Workflow with dynamic entry point detection
- **lib/EntryPoints.groovy** - Entry point registry and validation
- **DESIGN.md** - System design document
- **PLUGIN_PROPOSAL.md** - Proposal for Nextflow plugin

## Why It's Separate

This is an **advanced pattern** that adds complexity. Most users should start with the simple examples in the parent directory.

## Quick Start

```bash
cd advanced_entrypoints

# Show available entry points
nextflow run main.nf --help

# Auto-detect entry point
nextflow run main.nf -stub --outdir results

# Force specific entry point
nextflow run main.nf -stub --entry_point RUN_TOOLS --outdir results
```

## Features

✅ **Automatic entry point detection** - Detects which stage to run based on existing files
✅ **Dependency validation** - Ensures required files exist before starting
✅ **Clear error messages** - Shows exactly what's missing
✅ **Manual override** - Specify entry point explicitly with `--entry_point`

## Documentation

- [DESIGN.md](DESIGN.md) - How the system works
- [PLUGIN_PROPOSAL.md](PLUGIN_PROPOSAL.md) - Why this should be a Nextflow plugin

## Back to Basics

If this seems too complex, go back to the parent directory for simple examples:

```bash
cd ..
nextflow run main.nf -stub --outdir results
```

See [../INDEX.md](../INDEX.md) for all available examples.
