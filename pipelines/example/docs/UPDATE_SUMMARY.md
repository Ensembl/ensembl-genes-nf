# Structure Update Summary

The example pipeline has been reorganized for clarity.

## What Changed

### Before
- `main.nf` - Complex entry point system
- All files in root directory

### After
- **`main.nf`** - Simple workflow (two subworkflows)
- **`advanced_entrypoints/`** - Entry point system moved here

## New Structure

```
pipelines/example/
├── main.nf                          ← SIMPLE: Just runs two subworkflows
├── workflows/                       ← Learning examples
│   ├── simple_workflow.nf          
│   ├── minimal_example.nf
│   └── two_subworkflows.nf
├── subworkflows/                    ← Reusable components
│   ├── run_tools.nf
│   ├── combine_and_count.nf
│   └── ...
├── modules/                         ← Individual processes
│   └── ...
└── advanced_entrypoints/            ← ADVANCED: Entry point system
    ├── main.nf                      ← Dynamic entry points
    ├── lib/EntryPoints.groovy
    ├── DESIGN.md
    └── PLUGIN_PROPOSAL.md
```

## Quick Start (Updated)

### Simple Workflow (Recommended Start)
```bash
# Run the main workflow - simple and straightforward
nextflow run main.nf -stub --outdir results
```

### Advanced Entry Points
```bash
# Only if you need dynamic entry points
cd advanced_entrypoints
nextflow run main.nf --help
nextflow run main.nf -stub --outdir results
```

## Why This Change?

**Better learning experience:**
- New users see simple examples first
- Advanced patterns are clearly separated
- Less overwhelming for beginners

**Clearer organization:**
- Simple workflows in root
- Advanced features in subdirectory
- Documentation updated accordingly

## Migration

No breaking changes - all files still work, just in new locations:

- Entry point system: `main.nf` → `advanced_entrypoints/main.nf`
- Entry point lib: `lib/` → `advanced_entrypoints/lib/`
- Design docs moved to `advanced_entrypoints/`

## Documentation Updates

All documentation has been updated:
- [README.md](README.md) - Updated to describe simple main.nf
- [QUICK_START.md](QUICK_START.md) - Updated paths
- [INDEX.md](INDEX.md) - Updated structure
- [advanced_entrypoints/README.md](advanced_entrypoints/README.md) - New intro

