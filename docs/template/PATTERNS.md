# Common Nextflow Subworkflow Patterns

This guide shows common patterns you'll use when building Nextflow pipelines with subworkflows.

## Pattern 1: Parallel Execution

**Use case:** Run multiple processes on the same input in parallel.

```groovy
workflow PARALLEL_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // All processes run in parallel on the same input
    PROCESS_A(input_ch)
    PROCESS_B(input_ch)
    PROCESS_C(input_ch)

    // Combine all outputs
    all_results = PROCESS_A.out.results
        .concat(PROCESS_B.out.results)
        .concat(PROCESS_C.out.results)

    emit:
    results = all_results
}
```

**Example:** Running multiple QC tools on the same sequencing data.

---

## Pattern 2: Sequential Pipeline

**Use case:** Chain processes where each depends on the previous output.

```groovy
workflow SEQUENTIAL_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(raw_data) ]

    main:
    // Step 1
    TRIM(input_ch)

    // Step 2 uses Step 1 output
    ALIGN(TRIM.out.trimmed)

    // Step 3 uses Step 2 output
    SORT(ALIGN.out.bam)

    emit:
    bam    = SORT.out.sorted
    stats  = TRIM.out.stats.concat(ALIGN.out.stats)
}
```

**Example:** Read trimming → alignment → sorting.

---

## Pattern 3: Branch and Merge

**Use case:** Split data for different processing, then combine results.

```groovy
workflow BRANCH_MERGE_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // Split based on metadata
    input_ch.branch { meta, file ->
        longread: meta.type == 'long_read'
        shortread: meta.type == 'short_read'
    }
    .set { branched }

    // Different processing per branch
    LONGREAD_PROCESS(branched.longread)
    SHORTREAD_PROCESS(branched.shortread)

    // Merge results
    all_results = LONGREAD_PROCESS.out.results
        .concat(SHORTREAD_PROCESS.out.results)

    emit:
    results = all_results
}
```

**Example:** Different analysis for DNA vs RNA samples.

---

## Pattern 4: Conditional Execution

**Use case:** Run process only if condition is met.

```groovy
workflow CONDITIONAL_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // Always run
    PROCESS_A(input_ch)

    // Only run if meta.run_optional == true
    optional_input = PROCESS_A.out.results
        .filter { meta, file -> meta.run_optional }

    OPTIONAL_PROCESS(optional_input)

    // Combine outputs (handling empty channels)
    all_results = PROCESS_A.out.results
        .mix(OPTIONAL_PROCESS.out.results.ifEmpty([]))

    emit:
    results = all_results
}
```

**Example:** Optional quality filtering step.

---

## Pattern 5: Grouping and Aggregation

**Use case:** Collect multiple files per sample, then process together.

```groovy
workflow GROUP_AGGREGATE_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]
             // Multiple files per meta.id

    main:
    // Group all files by sample ID
    grouped = input_ch.groupTuple()
    // Now: [ val(meta), [file1, file2, file3, ...] ]

    // Process all files together
    AGGREGATE(grouped)

    emit:
    aggregated = AGGREGATE.out.results
}
```

**Example:** Combining multiple lanes of sequencing data.

---

## Pattern 6: Map and Collect

**Use case:** Process each item individually, then collect all results.

```groovy
workflow MAP_COLLECT_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // Process each file individually
    INDIVIDUAL_PROCESS(input_ch)

    // Collect all results into one channel
    all_files = INDIVIDUAL_PROCESS.out.results
        .collect()
    // Now: [ [file1, file2, file3, ...] ]

    // Process the collection
    SUMMARY(all_files)

    emit:
    individual = INDIVIDUAL_PROCESS.out.results
    summary    = SUMMARY.out.report
}
```

**Example:** Generate per-sample reports, then create overall summary.

---

## Pattern 7: Paired Inputs

**Use case:** Process paired files (e.g., R1/R2 reads) together.

```groovy
workflow PAIRED_PATTERN {
    take:
    reads_ch  // channel: [ val(meta), [path(R1), path(R2)] ]

    main:
    // Process paired files
    PAIRED_PROCESS(reads_ch)

    emit:
    results = PAIRED_PROCESS.out.results
}

// Creating paired channel from separate R1/R2 files:
Channel
    .fromFilePairs("data/*_R{1,2}.fastq.gz")
    .map { sample_id, files ->
        def meta = [id: sample_id]
        [meta, files]
    }
    .set { reads_ch }
```

**Example:** Processing paired-end sequencing reads.

---

## Pattern 8: Join Two Channels

**Use case:** Combine data from two different sources by sample ID.

```groovy
workflow JOIN_PATTERN {
    take:
    bam_ch       // channel: [ val(meta), path(bam) ]
    reference_ch // channel: [ val(meta), path(ref) ]

    main:
    // Join by meta.id
    joined = bam_ch.join(reference_ch)
    // Now: [ val(meta), path(bam), path(ref) ]

    PROCESS_WITH_REF(joined)

    emit:
    results = PROCESS_WITH_REF.out.results
}
```

**Example:** Joining sample data with sample-specific reference.

---

## Pattern 9: Mix Multiple Inputs

**Use case:** Combine outputs from different processes into one channel.

```groovy
workflow MIX_PATTERN {
    take:
    input_a  // channel: [ val(meta), path(file) ]
    input_b  // channel: [ val(meta), path(file) ]

    main:
    PROCESS_A(input_a)
    PROCESS_B(input_b)

    // Mix both outputs into one channel
    mixed = PROCESS_A.out.results
        .mix(PROCESS_B.out.results)

    DOWNSTREAM(mixed)

    emit:
    results = DOWNSTREAM.out.results
}
```

**Example:** Combining samples from different sequencing runs.

---

## Pattern 10: Multi-Output Process

**Use case:** Process produces multiple outputs that go to different downstream processes.

```groovy
workflow MULTI_OUTPUT_PATTERN {
    take:
    input_ch  // channel: [ val(meta), path(file) ]

    main:
    // Process with multiple outputs
    MULTI_OUTPUT_PROCESS(input_ch)

    // Different outputs to different processes
    PROCESS_A(MULTI_OUTPUT_PROCESS.out.output_a)
    PROCESS_B(MULTI_OUTPUT_PROCESS.out.output_b)

    emit:
    results_a = PROCESS_A.out.results
    results_b = PROCESS_B.out.results
}
```

**Example:** Alignment produces BAM file and unmapped reads separately.

---

## Pattern 11: Error Handling

**Use case:** Continue pipeline even if some samples fail.

```groovy
// In process definition
process PROCESS_WITH_ERRORS {
    errorStrategy 'ignore'  // or 'retry', 'finish'
    maxRetries 3

    // ... rest of process
}

workflow ERROR_HANDLING_PATTERN {
    take:
    input_ch

    main:
    PROCESS_WITH_ERRORS(input_ch)

    // Filter successful outputs
    successful = PROCESS_WITH_ERRORS.out.results
        .filter { meta, file -> file.exists() }

    emit:
    results = successful
}
```

---

## Pattern 12: Cross Product (Cartesian)

**Use case:** Run every combination of two sets.

```groovy
workflow CROSS_PATTERN {
    take:
    samples_ch   // channel: [ val(meta), path(sample) ]
    params_ch    // channel: [ val(param_set) ]

    main:
    // Create all combinations
    combinations = samples_ch.combine(params_ch)
    // Now: [ val(meta), path(sample), val(param_set) ]

    PROCESS_COMBO(combinations)

    emit:
    results = PROCESS_COMBO.out.results
}
```

**Example:** Testing multiple parameter sets on each sample.

---

## Best Practices

### 1. **Use meaningful names**
```groovy
// Good
aligned_bam = ALIGN.out.bam
sorted_bam = SORT.out.sorted

// Bad
output1 = ALIGN.out.bam
x = SORT.out.sorted
```

### 2. **Document channel structure**
```groovy
workflow MY_WORKFLOW {
    take:
    input_ch  // channel: [ val(meta), [path(R1), path(R2)] ]
              //   meta: [id: sample_id, single_end: false]

    emit:
    results  // channel: [ val(meta), path(bam) ]
}
```

### 3. **Handle empty channels**
```groovy
// Prevent errors from empty channels
optional_results = OPTIONAL.out.results.ifEmpty([])
all_results = required_results.mix(optional_results)
```

### 4. **Use consistent meta structure**
```groovy
// Standard meta format
meta = [
    id: 'sample_name',      // Required: sample identifier
    single_end: false,      // Optional: sequencing type
    strandedness: 'forward' // Optional: other metadata
]
```

### 5. **Test subworkflows independently**
```groovy
// Create minimal test workflow
include { MY_SUBWORKFLOW } from './subworkflows/my_subworkflow'

workflow {
    test_input = Channel.of(
        [[id: 'test1'], file('test1.txt')],
        [[id: 'test2'], file('test2.txt')]
    )

    MY_SUBWORKFLOW(test_input)
}
```

---

## Common Mistakes to Avoid

- **Forgetting to emit outputs**
```groovy
// Missing emit block - outputs not accessible!
workflow BROKEN {
    take: input_ch
    main:
    PROCESS(input_ch)
    // Missing: emit: results = PROCESS.out.results
}
```

- **Modifying meta in place**
```groovy
// Bad: modifies original meta
meta.new_field = 'value'

// Good: create new meta
def new_meta = meta + [new_field: 'value']
```

- **Not handling empty channels**
```groovy
// Bad: fails if optional_results is empty
all = required.concat(optional_results)

// Good: handle empty case
all = required.concat(optional_results.ifEmpty([]))
```

- **Hardcoding values**
```groovy
// Bad
publishDir "/home/user/results"

// Good
publishDir "${params.outdir}/results"
```

---

## See Also

- [Simple examples](README.md) - Basic workflow structure
- [Module template](../modules/minimal_example.nf) - Creating new processes

