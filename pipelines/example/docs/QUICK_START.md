# Quick Start Guide - Example Pipeline

Three levels of examples from simplest to most advanced.

## Level 1: Basic Workflows (Start Here!)

**Goal:** Understand basic Nextflow workflow structure

### Run the Simplest Example
```bash
cd pipelines/example
nextflow run workflows/simple_workflow.nf -stub --outdir simple_results
```

**What it does:** Runs one process on two samples

**Files to read:**
1. [workflows/simple_workflow.nf](workflows/simple_workflow.nf) - Main workflow
2. [modules/tool_a.nf](modules/tool_a.nf) - The process

**Next:** Try [workflows/minimal_example.nf](workflows/minimal_example.nf) to see subworkflows

---

## Level 2: Subworkflows

**Goal:** Learn to organize code with subworkflows

### Run Subworkflow Examples
```bash
# Two subworkflows chained together
nextflow run workflows/two_subworkflows.nf -stub --outdir results

# Minimal example
nextflow run workflows/minimal_example.nf -stub --outdir results
```

**What it does:**
- First subworkflow: Runs 4 tools in parallel
- Second subworkflow: Combines results and counts lines

**Files to read:**
1. [workflows/two_subworkflows.nf](workflows/two_subworkflows.nf) - Main workflow
2. [subworkflows/minimal_subworkflow_example.nf](subworkflows/minimal_subworkflow_example.nf) - Parallel pattern
3. [subworkflows/simple_sequential_example.nf](subworkflows/simple_sequential_example.nf) - Sequential pattern

**Guides:**
- [Workflow patterns guide](workflows/PATTERNS.md) - Common patterns
- [Workflows README](workflows/README.md) - Detailed docs

---

## Level 3: Dynamic Entry Points (Advanced)

**Goal:** Understand automatic workflow resumption

### Run with Entry Points
```bash
cd advanced_entrypoints

# First run: Full workflow
nextflow run main.nf -stub --outdir step1

# Second run: Automatic resume
rm -rf step1/combined step1/line_counts
nextflow run main.nf -stub --outdir step1
# ✓ Automatically detects existing tool outputs
# ✓ Skips to combine step
```

### Force Specific Entry Point
```bash
cd advanced_entrypoints

# Run only tools
nextflow run main.nf -stub --entry_point RUN_TOOLS --outdir results

# Run only combine (requires tool outputs)
nextflow run main.nf -stub --entry_point COMBINE_AND_COUNT --outdir results
```

### Get Help
```bash
cd advanced_entrypoints
nextflow run main.nf --help
```

**Files to read:**
1. [advanced_entrypoints/main.nf](advanced_entrypoints/main.nf) - Dynamic entry point workflow
2. [advanced_entrypoints/lib/EntryPoints.groovy](advanced_entrypoints/lib/EntryPoints.groovy) - Entry point registry
3. [advanced_entrypoints/DESIGN.md](advanced_entrypoints/DESIGN.md) - System design
4. [advanced_entrypoints/PLUGIN_PROPOSAL.md](advanced_entrypoints/PLUGIN_PROPOSAL.md) - Why this should be a plugin

---

## Quick Reference

### File Structure
```
pipelines/example/
├── workflows/              # Simple examples (START HERE)
│   ├── simple_workflow.nf
│   ├── minimal_example.nf
│   └── two_subworkflows.nf
├── subworkflows/           # Reusable subworkflows
│   ├── run_tools.nf
│   └── combine_and_count.nf
├── modules/                # Individual processes
│   ├── tool_a.nf
│   ├── tool_b.nf
│   └── ...
├── advanced_entrypoints/   # Advanced: Entry point system
│   ├── main.nf            # Dynamic entry points
│   └── lib/EntryPoints.groovy
└── main.nf                 # Simple: Two subworkflows
```

### Common Commands

```bash
# Run in stub mode (fast, for testing)
nextflow run <workflow.nf> -stub --outdir <dir>

# Run for real
nextflow run <workflow.nf> --outdir <dir>

# Generate workflow diagram
nextflow run <workflow.nf> -with-dag flowchart.html

# Resume from cached results
nextflow run <workflow.nf> -resume

# See workflow report
nextflow run <workflow.nf> -with-report report.html
```

### Learning Path

1. **Day 1:** Run `simple_workflow.nf` and understand the code
2. **Day 2:** Run `minimal_example.nf` and see how subworkflows connect
3. **Day 3:** Read [PATTERNS.md](workflows/PATTERNS.md) and try examples
4. **Day 4:** Explore the entry point system in `main.nf`
5. **Day 5:** Read [DESIGN.md](DESIGN.md) to understand advanced concepts

---

## Creating Your Own

### 1. Create a New Process

Copy the template:
```bash
cp modules/minimal_example.nf modules/my_tool.nf
```

Edit and replace TODOs:
- Process name
- Container URL
- Input/output channels
- Command to run

### 2. Create a New Subworkflow

```groovy
include { MY_TOOL } from '../modules/my_tool'

workflow MY_SUBWORKFLOW {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    MY_TOOL(input_ch)

    emit:
    results  = MY_TOOL.out.results
    versions = MY_TOOL.out.versions
}
```

### 3. Create a New Workflow

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MY_SUBWORKFLOW } from '../subworkflows/my_subworkflow'

params.outdir = 'results'

workflow {
    // Create input
    Channel.of([id: 'sample1'])
        .map { meta -> [meta, file('input.txt')] }
        .set { input_ch }

    // Run subworkflow
    MY_SUBWORKFLOW(input_ch)
}
```

---

## Getting Help

- **Workflow basics:** [workflows/README.md](workflows/README.md)
- **Common patterns:** [workflows/PATTERNS.md](workflows/PATTERNS.md)
- **Entry point system:** [DESIGN.md](DESIGN.md)
- **Plugin proposal:** [PLUGIN_PROPOSAL.md](PLUGIN_PROPOSAL.md)
- **Nextflow docs:** https://nextflow.io/docs/latest/

---

## Tips

✅ Always start in `-stub` mode for fast iteration
✅ Use descriptive channel names (not `ch1`, `ch2`)
✅ Document your channel structure with comments
✅ Test each subworkflow independently
✅ Use `groupTuple()` to collect files per sample

❌ Don't hardcode paths (use `params`)
❌ Don't modify `meta` in place (create new)
❌ Don't forget `emit:` block in subworkflows
❌ Don't forget to handle empty channels

---

**Ready to dive in?** Start with [workflows/simple_workflow.nf](workflows/simple_workflow.nf)!
