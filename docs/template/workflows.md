# Workflow Design Guide

How to design and structure workflows and subworkflows in this repository.

## Workflow Philosophy

**A useful default is: workflows orchestrate and modules execute.** This keeps
the dataflow easy to follow, but small pipelines may reasonably combine a few
layers when introducing more files would add noise.

## Keep parameter ownership at the entrypoint

The entrypoint (`main.nf`) is usually the best place to own named pipeline
parameters. It validates
the schema, resolves parameter-dependent files and defaults, and translates
those values into explicit channels or values before calling the workflow
layer. This is the boundary between the user's configuration and the
pipeline's dataflow.

For reusable named workflows and subworkflows, it is usually clearer not to
reach into `params` for inputs. Instead, declare dependencies in `take:` and
pass them as arguments. For example:

```groovy
// main.nf
workflow {
    validateParameters()
    samples_ch = channel.fromPath(params.input)
    ANALYSIS(samples_ch, params.mode, file(params.reference))
}

// subworkflows/analysis.nf
workflow ANALYSIS {
    take:
    samples_ch
    mode
    reference

    main:
    if (mode == 'combined') {
        TOOL_A(samples_ch, reference)
        TOOL_B(samples_ch, reference)
    }
}
```

This reduces coupling between a reusable subworkflow and a particular
parameter name, schema, or command-line interface. Another pipeline can reuse
`ANALYSIS` with a different manifest or parameter schema, provided it supplies
the same explicit inputs. It also makes dependencies visible in the call site,
which can make testing, documentation, and later changes easier. Direct
`params` access can still be reasonable for genuinely pipeline-wide settings;
the important distinction is whether the dependency is intentional and clear.

The intended direction is therefore:

```text
named parameters → main.nf → workflows → subworkflows → modules
```

Each layer can narrow configuration into explicit dataflow. Parameters are
often easiest to manage at the edge, while reusable components generally work
best with declared inputs.

### Key Principles

1. **Workflows coordinate** - They connect modules and manage data flow
2. **Keep responsibilities visible** - Put logic where it is easiest to test and understand
3. **Make subworkflows reusable when useful** - Reuse is a benefit, not a requirement
4. **Use a clear hierarchy when the pipeline needs it** - Main workflow → subworkflows → modules

## Workflow vs Subworkflow

### When to Use a Workflow

A **workflow** is your main entry point:

```groovy
// main.nf or my_pipeline.nf
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SUBWORKFLOW_A } from './subworkflows/subworkflow_a'
include { SUBWORKFLOW_B } from './subworkflows/subworkflow_b'

params.outdir = 'results'

workflow {
    // Create input
    input_ch = channel.fromPath(params.input)

    // Orchestrate subworkflows
    SUBWORKFLOW_A(input_ch)
    SUBWORKFLOW_B(SUBWORKFLOW_A.out.results)
}
```

**Typical characteristics**:
- Entry point for `nextflow run`
- No `take` or `emit` blocks (unnamed workflow)
- Defines pipeline-level parameters
- Coordinates high-level workflow logic

### When to Use a Subworkflow

A **subworkflow** is a reusable component:

```groovy
// subworkflows/process_samples.nf
include { TOOL_A } from '../modules/tool_a'
include { TOOL_B } from '../modules/tool_b'

workflow PROCESS_SAMPLES {
    take:
    samples_ch  // channel: [ val(meta), path(file) ]

    main:
    TOOL_A(samples_ch)
    TOOL_B(TOOL_A.out.processed)

    emit:
    results  = TOOL_B.out.results
    stats    = TOOL_A.out.stats.concat(TOOL_B.out.stats)
    versions = TOOL_A.out.versions.concat(TOOL_B.out.versions)
}
```

**Typical characteristics**:
- Named workflow with `take` and `emit` blocks
- Reusable across multiple pipelines
- Encapsulates a logical unit of work
- Clear inputs and outputs

**A subworkflow is especially useful when**:
- Logic is reused across multiple pipelines
- Grouping related processes makes sense conceptually
- You want to test a component independently
- It represents a distinct analysis stage (QC, alignment, quantification)

## Subworkflow Design Patterns

### Pattern 1: Parallel Execution

**Use case**: Run multiple independent processes on the same input

```groovy
workflow PARALLEL_TOOLS {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // All processes receive same input, run in parallel
    TOOL_A(input_ch)
    TOOL_B(input_ch)
    TOOL_C(input_ch)

    // Collect outputs
    all_results = TOOL_A.out.results
        .concat(TOOL_B.out.results)
        .concat(TOOL_C.out.results)

    emit:
    results  = all_results
    versions = TOOL_A.out.versions
        .concat(TOOL_B.out.versions)
        .concat(TOOL_C.out.versions)
}
```

**Why this pattern**:
- Maximum parallelization - all tools run simultaneously
- Independent results - each tool produces separate output
- Common in QC (multiple QC tools) or multi-tool analysis

### Pattern 2: Sequential Pipeline

**Use case**: Chain processes where each depends on previous output

```groovy
workflow SEQUENTIAL_ANALYSIS {
    take:
    raw_data  // channel: [ val(meta), path(reads) ]

    main:
    // Each step uses output of previous
    TRIM(raw_data)
    ALIGN(TRIM.out.trimmed)
    SORT(ALIGN.out.bam)
    INDEX(SORT.out.sorted)

    emit:
    bam      = SORT.out.sorted
    bai      = INDEX.out.index
    versions = TRIM.out.versions
        .concat(ALIGN.out.versions)
        .concat(SORT.out.versions)
        .concat(INDEX.out.versions)
}
```

**Why this pattern**:
- Enforces dependencies - steps run in correct order
- Clear data flow - easy to understand pipeline logic
- Common for preprocessing pipelines

### Pattern 3: Conditional Branching

**Use case**: Different processing based on data characteristics

```groovy
workflow CONDITIONAL_PROCESSING {
    take:
    samples_ch  // channel: [ val(meta), path(file) ]

    main:
    // Branch by sample type
    samples_ch
        .branch { meta, file ->
            long_read:  meta.read_type == 'long'
            short_read: meta.read_type == 'short'
        }
        .set { branched }

    // Different processing per branch
    LONGREAD_ALIGN(branched.long_read)
    SHORTREAD_ALIGN(branched.short_read)

    // Merge results back
    all_aligned = LONGREAD_ALIGN.out.bam
        .mix(SHORTREAD_ALIGN.out.bam)

    emit:
    aligned = all_aligned
    versions = LONGREAD_ALIGN.out.versions
        .mix(SHORTREAD_ALIGN.out.versions)
}
```

**Why this pattern**:
- Handles heterogeneous data appropriately
- Each branch uses specialized tools
- Results reconverge for downstream processing

### Pattern 4: Aggregation

**Use case**: Collect per-sample results, process together

```groovy
workflow AGGREGATE_ANALYSIS {
    take:
    per_sample_results  // channel: [ val(meta), path(file) ]

    main:
    // Group all files by sample
    grouped = per_sample_results.groupTuple()

    // Process per sample
    PER_SAMPLE_SUMMARY(grouped)

    // Collect all samples for cohort analysis
    all_samples = PER_SAMPLE_SUMMARY.out.summary.collect()

    // Cohort-level analysis
    COHORT_ANALYSIS(all_samples)

    emit:
    per_sample = PER_SAMPLE_SUMMARY.out.summary
    cohort     = COHORT_ANALYSIS.out.results
    versions   = PER_SAMPLE_SUMMARY.out.versions
        .concat(COHORT_ANALYSIS.out.versions)
}
```

**Why this pattern**:
- `groupTuple()` - Collects files with same meta.id
- `collect()` - Waits for all samples, creates single channel
- Common for multi-sample comparisons or cohort statistics

## Workflow Structure Best Practices

### Input Channel Creation

**Create channels at the top of workflow**:

```groovy
workflow {
    // Good: Clear input creation
    Channel
        .fromPath(params.input_files)
        .map { file ->
            def meta = [id: file.baseName]
            [meta, file]
        }
        .set { input_ch }

    // Then use in subworkflows
    MY_SUBWORKFLOW(input_ch)
}
```

**Why**:
- All input creation in one place
- Easy to modify input handling
- Clear what data enters the pipeline

### Subworkflow Outputs

**Emit what downstream needs**:

```groovy
workflow MY_SUBWORKFLOW {
    take:
    input_ch

    main:
    PROCESS_A(input_ch)
    PROCESS_B(PROCESS_A.out.results)

    emit:
    // Primary outputs
    results = PROCESS_B.out.results

    // Intermediate outputs (if needed downstream)
    intermediate = PROCESS_A.out.results

    // Emit versions for all processes when version tracking is part of the pipeline contract
    versions = PROCESS_A.out.versions
        .concat(PROCESS_B.out.versions)
}
```

**Why**:
- Named emits are self-documenting
- Flexible for different use cases
- Version tracking is standard

### Channel Documentation

**Document channel structure with comments**:

```groovy
workflow ALIGNMENT {
    take:
    reads_ch  // channel: [ val(meta), [ path(R1), path(R2) ] ]
              //   meta: [id: sample_id, single_end: false]
    ref       // path: reference genome fasta

    main:
    // Process implementation

    emit:
    bam       // channel: [ val(meta), path(bam) ]
    bai       // channel: [ val(meta), path(bai) ]
    stats     // channel: [ val(meta), path(stats_txt) ]
    versions  // path: versions.yml
}
```

**Why**:
- Clear expectations for inputs
- Easy to understand data flow
- Helps with debugging

## Common Workflow Patterns

### Pattern: Combining Multiple Outputs

```groovy
main:
TOOL_A(input_ch)
TOOL_B(input_ch)
TOOL_C(input_ch)

// Combine outputs from parallel processes
combined = TOOL_A.out.results
    .concat(TOOL_B.out.results)
    .concat(TOOL_C.out.results)
```

**When to use**: Parallel processes producing similar outputs

### Pattern: Joining Data from Different Sources

```groovy
main:
PROCESS_A(input_a)  // Produces: [ meta, file_a ]
PROCESS_B(input_b)  // Produces: [ meta, file_b ]

// Join by meta.id
joined = PROCESS_A.out.results
    .join(PROCESS_B.out.results)  // Produces: [ meta, file_a, file_b ]

DOWNSTREAM(joined)
```

**When to use**: Need to combine files from different processes by sample

### Pattern: Collecting All Samples

```groovy
main:
PER_SAMPLE(input_ch)

// Wait for all samples, create single channel
all_samples = PER_SAMPLE.out.results
    .collect()  // Produces: [ [file1, file2, file3, ...] ]

MULTI_SAMPLE_ANALYSIS(all_samples)
```

**When to use**: Cohort analysis, multi-sample comparisons

### Pattern: Optional Processes

```groovy
main:
REQUIRED_STEP(input_ch)

// Filter based on parameter
optional_input = REQUIRED_STEP.out.results
    .filter { meta, file -> !params.skip_optional }

OPTIONAL_STEP(optional_input)

// Combine with ifEmpty to handle missing data
all_results = REQUIRED_STEP.out.results
    .mix(OPTIONAL_STEP.out.results.ifEmpty([]))
```

**When to use**: Optional QC, conditional analysis steps

## Workflow Organization

### File Structure

```
pipeline_name/
├── main.nf                    # Entry point workflow
├── nextflow.config            # Pipeline configuration
├── workflows/                 # Alternative entry points
│   ├── workflow_variant_a.nf
│   └── workflow_variant_b.nf
├── subworkflows/              # Reusable components
│   ├── qc.nf
│   ├── alignment.nf
│   └── quantification.nf
└── modules/                   # Individual processes
    ├── fastqc.nf
    ├── trim_galore.nf
    └── star.nf
```

**Hierarchy principle**: main.nf → subworkflows → modules

### Naming Conventions

**Workflows**:
- Entry point: `main.nf` or descriptive name `rnaseq.nf`
- Alternatives: `workflow_<variant>.nf`

**Subworkflows**:
- Uppercase names: `ALIGNMENT`, `QC`, `QUANTIFICATION`
- Verb or noun describing function
- File names: lowercase with underscores `alignment.nf`

**Include statements**:
```groovy
// Subworkflows (uppercase in code)
include { ALIGNMENT } from './subworkflows/alignment'
include { QC } from './subworkflows/qc'

// Modules (uppercase in code)
include { FASTQC } from './modules/fastqc'
include { TRIMGALORE } from './modules/trim_galore'
```

## Testing Workflows

### Stub Mode Testing

**Use stub mode to validate structure**:

```bash
# Test workflow logic without running tools
nextflow run main.nf -stub --outdir test_results
```

**Why**:
- Fast - seconds instead of hours
- Validates channel connections
- Tests branching logic
- Catches naming issues

### Workflow Visualization

**Generate DAG to understand structure**:

```bash
# Create visual workflow diagram
nextflow run main.nf -with-dag workflow.html

# Or as SVG
nextflow run main.nf -with-dag workflow.svg
```

**Why**:
- See parallelization opportunities
- Verify dependencies
- Documentation for team

## Anti-Patterns to Avoid

### Don't: Put Logic in Workflows

```groovy
// Bad: Business logic in workflow
workflow {
    input_ch.map { meta, file ->
        if (meta.type == 'A') {
            // Complex processing logic here
        }
    }
}

// Good: Logic in module, workflow coordinates
workflow {
    PROCESS_SAMPLES(input_ch)
}
```

### Don't: Hardcode Values

```groovy
// Bad: Hardcoded paths and values
workflow {
    reference = file("/hardcoded/path/to/ref.fa")
    ALIGN(input_ch, reference)
}

// Good: Use parameters
workflow {
    reference = file(params.reference)
    ALIGN(input_ch, reference)
}
```

### Don't: Create Monolithic Subworkflows

```groovy
// Bad: Everything in one subworkflow
workflow DO_EVERYTHING {
    // QC + alignment + quantification + statistics
    // 30 processes all in one subworkflow
}

// Good: Logical separation
workflow QC { /* ... */ }
workflow ALIGNMENT { /* ... */ }
workflow QUANTIFICATION { /* ... */ }
```

### Don't: Forget Error Handling

```groovy
// Bad: No handling of empty channels
workflow {
    OPTIONAL_STEP(input_ch)
    DOWNSTREAM(OPTIONAL_STEP.out.results)  // Fails if empty!
}

// Good: Handle empty channels
workflow {
    OPTIONAL_STEP(input_ch)
    results = OPTIONAL_STEP.out.results.ifEmpty([])
    DOWNSTREAM(results)
}
```

### Don't: Duplicate Subworkflows

```groovy
// Bad: Copy-paste similar subworkflows
workflow ALIGN_SAMPLE_A { /* ... */ }
workflow ALIGN_SAMPLE_B { /* ... */ }  // Almost identical

// Good: One parameterized subworkflow
workflow ALIGN {
    take:
    samples_ch

    main:
    // Handle both cases with branching or parameters
}
```

## Best Practices Summary

**Workflow Design**:
- Keep workflows thin - coordinate, don't implement
- Use subworkflows for reusable logic
- Document channel structures with comments
- Test with stub mode early and often

**Subworkflow Design**:
- Clear take/emit blocks with types
- One logical unit of work per subworkflow
- Emit versions when they are part of the pipeline's reproducibility contract
- Design for reusability

**Channel Operations**:
- Create inputs at top of workflow
- Use named emits for clarity
- Handle empty channels with `ifEmpty()`
- Join/combine thoughtfully

**Organization**:
- Follow hierarchy: workflow → subworkflow → module
- Consistent naming conventions
- Separate configuration from logic
- Document with comments and DAGs

## Further Reading

- [MODULES.md](MODULES.md) - Module design principles
- [PATTERNS.md](PATTERNS.md) - 12 common workflow patterns
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration strategy
- [Nextflow Workflow Docs](https://nextflow.io/docs/latest/workflow.html)
