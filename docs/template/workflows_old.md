# Simple Workflow Examples

This directory contains minimal examples showing basic Nextflow workflow and subworkflow patterns, **without** the complexity of the entry point system.

## Examples

### 1. Simple Workflow (`simple_workflow.nf`)

The absolute minimum: one workflow, one process.

```bash
nextflow run workflows/simple_workflow.nf -stub --outdir simple_results
```

**Structure:**
```
workflow
  └── TOOL_A process
```

**Learn:**
- Basic workflow structure
- Creating input channels
- Running a single process

---

### 2. Minimal Example (`minimal_example.nf`)

Basic example with two subworkflows chained together.

```bash
nextflow run workflows/minimal_example.nf -stub --outdir results
```

**Structure:**
```
workflow
  ├── RUN_TOOLS subworkflow
  │   ├── TOOL_A
  │   ├── TOOL_B
  │   ├── TOOL_C
  │   └── TOOL_D
  └── COMBINE_AND_COUNT subworkflow
      ├── COMBINE_OUTPUTS
      └── COUNT_LINES
```

**Learn:**
- Including subworkflows
- Passing channels between subworkflows
- Using subworkflow outputs

---

### 3. Two Subworkflows (`two_subworkflows.nf`)

Shows how to create and chain custom subworkflows.

```bash
nextflow run workflows/two_subworkflows.nf -stub --outdir results
```

**Structure:**
```
workflow
  ├── MINIMAL_SUBWORKFLOW_EXAMPLE
  │   ├── TOOL_A
  │   └── TOOL_B
  └── SIMPLE_SEQUENTIAL_EXAMPLE
      ├── COMBINE_OUTPUTS
      └── COUNT_LINES
```

**Learn:**
- Creating custom subworkflows
- Parallel process execution in subworkflows
- Sequential process execution in subworkflows

---

## Subworkflow Examples

### Minimal Subworkflow (`minimal_subworkflow_example.nf`)

Shows basic subworkflow structure with parallel execution.

**Pattern:**
```groovy
workflow MINIMAL_SUBWORKFLOW_EXAMPLE {
    take:
    input_channel

    main:
    TOOL_A(input_channel)
    TOOL_B(input_channel)

    all_results = TOOL_A.out.results.concat(TOOL_B.out.results)

    emit:
    results  = all_results
    versions = TOOL_A.out.versions.concat(TOOL_B.out.versions)
}
```

**Use when:** You need to run multiple processes in parallel on the same input.

---

### Sequential Subworkflow (`simple_sequential_example.nf`)

Shows sequential process execution (output of one → input to next).

**Pattern:**
```groovy
workflow SIMPLE_SEQUENTIAL_EXAMPLE {
    take:
    input_files

    main:
    COMBINE_OUTPUTS(input_files.groupTuple())
    COUNT_LINES(COMBINE_OUTPUTS.out.combined)

    emit:
    combined = COMBINE_OUTPUTS.out.combined
    counts   = COUNT_LINES.out.counts
    versions = COMBINE_OUTPUTS.out.versions.concat(COUNT_LINES.out.versions)
}
```

**Use when:** Processes must run in order (one depends on the other's output).

---

## Key Concepts

### 1. Workflow Structure

```groovy
workflow {
    // 1. Create input
    Channel.of(...).set { input_ch }

    // 2. Run subworkflow
    MY_SUBWORKFLOW(input_ch)

    // 3. Use outputs
    NEXT_STEP(MY_SUBWORKFLOW.out.results)
}
```

### 2. Subworkflow Structure

```groovy
workflow MY_SUBWORKFLOW {
    take:
    input_channel  // Define inputs

    main:
    // Run processes
    PROCESS_A(input_channel)

    emit:
    results = PROCESS_A.out.results  // Define outputs
}
```

### 3. Including Modules/Subworkflows

```groovy
// Include a single process
include { TOOL_A } from '../modules/tool_a'

// Include multiple processes
include { TOOL_A; TOOL_B } from '../modules/tools'

// Include a subworkflow
include { RUN_TOOLS } from '../subworkflows/run_tools'
```

### 4. Channel Operations

```groovy
// Combine outputs from multiple processes
all_results = TOOL_A.out.results
    .concat(TOOL_B.out.results)

// Group by sample ID (meta.id)
grouped = channel.groupTuple()

// Pass through unchanged
NEXT_PROCESS(PREVIOUS.out.results)
```

---

## Running the Examples

All examples support the same parameters:

```bash
# Run with stub mode (fast, for testing)
nextflow run workflows/<workflow>.nf -stub --outdir <output_dir>

# Run for real
nextflow run workflows/<workflow>.nf --outdir <output_dir>

# See the workflow diagram
nextflow run workflows/<workflow>.nf -with-dag flowchart.html
```

---

## Progression

**Start here:**
1. `simple_workflow.nf` - Understand basic workflow structure
2. `minimal_subworkflow_example.nf` - Learn subworkflow patterns
3. `minimal_example.nf` - See how to chain subworkflows

**Then explore:**
- `two_subworkflows.nf` - Create custom subworkflows
- `../main.nf` - Advanced: Dynamic entry point system

---

## File Locations

```
pipelines/example/
├── workflows/                    # Main workflow files
│   ├── simple_workflow.nf       # Simplest example
│   ├── minimal_example.nf       # Basic chained subworkflows
│   └── two_subworkflows.nf      # Custom subworkflows
├── subworkflows/                # Reusable subworkflows
│   ├── run_tools.nf            # Parallel execution example
│   ├── combine_and_count.nf    # Sequential execution example
│   ├── minimal_subworkflow_example.nf
│   └── simple_sequential_example.nf
└── modules/                     # Individual processes
    ├── tool_a.nf
    ├── tool_b.nf
    ├── tool_c.nf
    ├── tool_d.nf
    ├── combine_outputs.nf
    └── count_lines.nf
```

---

## Tips

- **Keep workflows simple** - Main workflow should just orchestrate subworkflows
- **Subworkflows should be reusable** - Don't hardcode values
- **Use clear naming** - Process and subworkflow names should be descriptive
- **Document inputs/outputs** - Use comments in `take:` and `emit:` blocks
- **Test in stub mode** - Fast iteration during development
