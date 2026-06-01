# Getting Started with the Example Pipeline

This guide walks you through the example pipeline, starting with the simplest concepts and progressively building to more advanced patterns.

## Your First Workflow

Start by running a single process workflow to understand the basic structure.
```bash
cd pipelines/example
nextflow run workflows/simple_workflow.nf -stub --outdir simple_results
```

**What it does:** Runs one process on two samples

**Key concepts:**
- Basic workflow structure
- Creating input channels
- Running a single process

**Explore the code:**
- [workflows/simple_workflow.nf](workflows/simple_workflow.nf) - Main workflow
- [modules/tool_a.nf](modules/tool_a.nf) - The process definition

---

## Working with Subworkflows

Once you're comfortable with basic workflows, learn how to organize code into reusable subworkflows.
```bash
# Two subworkflows chained together
nextflow run workflows/two_subworkflows.nf -stub --outdir results

# Example with multiple subworkflows
nextflow run workflows/minimal_example.nf -stub --outdir results
```

**What it does:**
- First subworkflow: Runs tools in parallel
- Second subworkflow: Combines results and counts lines

**Key concepts:**
- Including and chaining subworkflows
- Passing channels between subworkflows
- Parallel and sequential execution patterns

**Explore the code:**
- [workflows/two_subworkflows.nf](workflows/two_subworkflows.nf) - Chaining custom subworkflows
- [subworkflows/minimal_subworkflow_example.nf](subworkflows/minimal_subworkflow_example.nf) - Parallel execution
- [subworkflows/simple_sequential_example.nf](subworkflows/simple_sequential_example.nf) - Sequential execution

**Additional resources:**
- [PATTERNS.md](PATTERNS.md) - 12 common Nextflow patterns
- [workflows.md](workflows.md) - Detailed workflow documentation

---

## Dynamic Entry Points

For more complex pipelines, you may want automatic workflow resumption based on existing outputs.
```bash
cd advanced_entrypoints

# See available entry points
nextflow run main.nf --help

# Automatic detection: runs full workflow if no outputs exist
nextflow run main.nf -stub --outdir results

# Automatic resume: detects existing outputs and skips completed steps
rm -rf results/combined results/line_counts
nextflow run main.nf -stub --outdir results
# Detects existing tool outputs and resumes from combine step

# Manual override: force a specific entry point
nextflow run main.nf -stub --entry_point RUN_TOOLS --outdir results
nextflow run main.nf -stub --entry_point COMBINE_AND_COUNT --outdir results
```

**Key concepts:**
- Automatic entry point detection based on existing outputs
- Dependency validation before execution
- Manual entry point override when needed
- Clear error messages for missing dependencies


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

### Suggested Learning Path

Follow this progression to build your understanding:

1. Run and explore `simple_workflow.nf` to understand basic workflow structure
2. Run `minimal_example.nf` to see how subworkflows connect
3. Review [PATTERNS.md](PATTERNS.md) for common patterns you'll use
4. Explore the entry point system in `advanced_entrypoints/`
5. Read the design documentation to understand advanced concepts

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

## Additional Resources

- **Workflow patterns:** [PATTERNS.md](PATTERNS.md) - 12 common patterns
- **Workflow examples:** [workflows.md](workflows.md) - Detailed examples
- **Entry point design:** [advanced_entrypoints/docs/advanced_design.md](../pipelines/example/advanced_entrypoints/docs/advanced_design.md)
- **Nextflow docs:** https://nextflow.io/docs/latest/
- **nf-core guidelines:** https://nf-co.re/docs/guidelines

---

## Tips

- Always start in `-stub` mode for fast iteration
- Use descriptive channel names (not `ch1`, `ch2`)
- Document your channel structure with comments
- Test each subworkflow independently
- Use `groupTuple()` to collect files per sample

- Don't hardcode paths (use `params`)
- Don't modify `meta` in place (create new)
- Don't forget `emit:` block in subworkflows
- Don't forget to handle empty channels

---

**Ready to dive in?** Start with [workflows/simple_workflow.nf](workflows/simple_workflow.nf)!
