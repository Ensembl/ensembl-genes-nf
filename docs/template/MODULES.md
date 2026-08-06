# Module Design Guide

How to design and structure Nextflow process modules in this repository.

## Module Philosophy

**A module is a single, reusable process** that wraps one primary bioinformatics tool. Modules should be:

- **Generic** - Work across different pipelines and use cases
- **Configurable** - Parameters come from config, not hardcoded
- **Consistent** - Follow standard patterns for inputs, outputs, and metadata
- **Self-documenting** - Include versions, clear naming, and structure

## Core Module Structure

```groovy
process TOOL_NAME {
    label 'process_medium'              // Resource allocation
    container "..."                     // Container definition
    tag "${meta.id}"                    // Process identification
    errorStrategy 'ignore'              // Error handling
    publishDir "${params.outdir}/..."   // Output location

    input:
    tuple val(meta), path(input)        // Metadata + files

    output:
    tuple val(meta), path("*.txt"),     emit: results
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tool_command ...
    """

    stub:
    """
    touch output.txt
    """
}
```

## Essential Components

### 1. Meta Map Pattern

**Why**: Track sample identity and attributes through the pipeline

**Structure**:
```groovy
meta = [
    id: 'sample_name',      // Required: unique sample identifier
    single_end: false,      // Optional: sequencing type
    // Add pipeline-specific fields as needed
]
```

**Usage in modules**:
```groovy
input:
tuple val(meta), path(reads)    // Meta travels with the data

output:
tuple val(meta), path("*.bam")  // Meta preserved in output

script:
def prefix = meta.id            // Use meta.id for file naming
"""
tool --sample ${meta.id} --input ${reads} --output ${prefix}.bam
"""
```

**Why this matters**:
- **Identity tracking** - Know which sample produced which output
- **Attribute propagation** - Sample properties flow through pipeline
- **Grouping operations** - Group outputs by `meta.id` using `groupTuple()`
- **Conditional logic** - Handle samples differently based on metadata

**Best practices**:
```groovy
- Always include meta.id
- Use meta for sample attributes, not file properties
- Don't modify meta in place, create new: meta + [new_field: value]
- Keep meta structure consistent across modules
```

### 2. Container Definitions

**Why**: Reproducibility and portability

**Where to define**:
```groovy
// Option 1: In module (explicit, visible)
process BLAST {
    container ""
}

// Option 2: In config (centralised, easier to update)
// nextflow.config
process {
    withName: 'BLAST' {
        container = "https://depot.galaxyproject.org/singularity/blast:2.12.0"
    }
}
```

**Best practices**:
```groovy
- Use specific version tags, not 'latest'
- Prefer Singularity for HPC environments
- Use BioContainers when available (depot.galaxyproject.org)
- Document container source in comments
```

**Example with justification**:
```groovy
process MINIMAP2 {
    // Using BioContainers for reproducibility
    // Version 2.24 required for long-read alignment features
    container "https://depot.galaxyproject.org/singularity/minimap2:2.24--h7132678_1"
}
```

### 3. Version Tracking

**Why**: Reproducibility - know exactly which tool versions produced results

**Required output**:
```groovy
output:
path "versions.yml", emit: versions
```

**Implementation in script**:
```groovy
script:
"""
tool_command --input data.fq --output results.txt

# Capture tool version
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool_name: \$(tool_name --version | sed 's/tool_name v//g')
END_VERSIONS
"""
```

**Why YAML format**:
- Structured and parseable
- Easily merged across processes
- Standard format used by nf-core
- Can be automatically collected into summary report

**Multiple tools**:
```groovy
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(samtools --version | head -n1 | sed 's/samtools //g')
    bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //g')
END_VERSIONS
```

### 4. Stub Blocks

**Why**: Fast workflow testing without running expensive computations

**Purpose**:
- Test workflow logic and connections
- Validate channel operations
- Check file naming and structure
- Develop pipelines quickly

**Implementation**:
```groovy
stub:
def prefix = task.ext.prefix ?: meta.id
"""
# Create all expected output files (empty or minimal)
touch ${prefix}.bam
touch ${prefix}.bai
touch ${prefix}_stats.txt

# Versions file (use known versions)
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tool_name: 2.1.0
END_VERSIONS
"""
```

**Running stub mode**:
```bash
nextflow run pipeline.nf -stub --outdir results
```

**Best practices**:
```groovy
- Create ALL output files that script creates
- Use same file naming logic as script block
- Keep stub simple - just touch/echo
- Include versions.yml with hardcoded versions
```

**Why this matters**:
- Workflows with 100+ processes can be tested in seconds
- Catch channel/naming issues early
- Safe to run repeatedly during development
- No need for test data or containers

## Module Design Patterns

### Pattern 1: Configurable Arguments via ext.args

**Why**: Keep modules generic and reusable

```groovy
// In module
script:
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: meta.id
"""
tool ${args} --input data.fq --output ${prefix}.txt
"""
```

```groovy
// In config
process {
    withName: 'TOOL' {
        ext.args = '--quality 30 --threads 4'
        ext.prefix = { "${meta.id}_processed" }
    }
}
```

**Benefits**:
- Module code never needs editing for parameter changes
- Different pipelines can use same module differently
- Configuration separate from implementation

### Pattern 2: Conditional Execution via ext.when

**Why**: Control process execution from config, not code

```groovy
// In module
when:
task.ext.when == null || task.ext.when

// In config
process {
    withName: 'OPTIONAL_QC' {
        ext.when = { !params.skip_qc && meta.requires_qc }
    }
}
```

**Benefits**:
- Modules don't need to know about skip parameters
- Complex logic kept in config
- Easy to enable/disable steps per pipeline

### Pattern 3: Dynamic File Naming

**Why**: Avoid name collisions, support multiple samples

```groovy
script:
def prefix = task.ext.prefix ?: meta.id
def suffix = task.ext.suffix ?: 'processed'
"""
tool --input ${input} --output ${prefix}_${suffix}.txt
"""
```

**Best practices**:
```groovy
- Always use meta.id in output names
- Use task.ext.prefix for customisation
- Check for input/output name conflicts
- Use consistent naming patterns across modules
```

### Pattern 4: Multiple Outputs with Emit Names

**Why**: Clear, named outputs for workflow connections

```groovy
output:
tuple val(meta), path("*.bam"),     emit: alignment
tuple val(meta), path("*.bai"),     emit: index
tuple val(meta), path("*_stats.txt"), emit: stats
path "versions.yml",                emit: versions
```

**Usage in workflow**:
```groovy
ALIGNER(input_ch)

// Use specific outputs
DOWNSTREAM_TOOL(ALIGNER.out.alignment)
STATISTICS(ALIGNER.out.stats)
```

**Benefits**:
- Self-documenting workflow connections
- Type safety (wrong output name = error)
- Clear data flow

### Pattern 5: Optional Outputs

**Why**: Handle tools that conditionally produce files

```groovy
output:
tuple val(meta), path("*.bam"),           emit: bam
tuple val(meta), path("*_report.html"),   emit: report, optional: true
tuple val(meta), path("*_warnings.txt"),  emit: warnings, optional: true
```

**When to use**:
- Reports only generated if issues found
- Outputs depend on input data characteristics
- Tool has optional features

## Common Patterns

### Paired-End vs Single-End Handling

```groovy
script:
def args = task.ext.args ?: ''
def prefix = meta.id
def is_paired = meta.single_end ? false : true
def input_files = is_paired ? "${reads[0]} ${reads[1]}" : "${reads}"
"""
aligner ${args} --input ${input_files} --output ${prefix}.bam
"""
```

### Reference Files (Shared Across Samples)

```groovy
input:
tuple val(meta), path(reads)    // Per-sample files
path reference                  // Shared reference (no meta)
path index                      // Shared index (no meta)
```

**Why separate**: Reference files shared by all samples, don't need meta tracking

### Multi-Command Processes

```groovy
// In config
process {
    withName: 'ALIGN_AND_SORT' {
        ext.args = '--very-sensitive'        // For aligner
        ext.args2 = '-@ 4'                   // For sorter
    }
}

// In module
script:
def args = task.ext.args ?: ''
def args2 = task.ext.args2 ?: ''
"""
aligner ${args} input.fq | sorter ${args2} > output.bam
"""
```

## Anti-Patterns to Avoid

- **Hardcoded parameters**
```groovy
// Bad
"""
tool --quality 30 --threads 4 input.fq
"""

// Good
"""
tool ${args} --threads ${task.cpus} input.fq
"""
```

- **Modifying meta in place**
```groovy
// Bad
meta.processed = true

// Good
def new_meta = meta + [processed: true]
```

- **Multiple primary tools in one module**
```groovy
// Bad - creates monolithic, inflexible modules
"""
align input.fq > aligned.bam
sort aligned.bam > sorted.bam
index sorted.bam
"""

// Good - separate modules, composable
// Use subworkflows to chain them
```

- **Forgetting stub block**
```groovy
// Bad - can't test workflow without running
script:
"""
expensive_computation
"""

// Good - rapid testing possible
script:
"""
expensive_computation
"""

stub:
"""
touch output.txt
"""
```

- **Missing versions.yml**
```groovy
// Bad - no reproducibility tracking
output:
tuple val(meta), path("*.txt"), emit: results

// Good - versions tracked
output:
tuple val(meta), path("*.txt"), emit: results
path "versions.yml",            emit: versions
```

## Module Template

See [modules/minimal_example.nf](../modules/minimal_example.nf) for a TODO-annotated template and [modules/reference_example.nf](../modules/reference_example.nf) for a comprehensive example with extensive comments.

## Best Practices Summary

- **Always include meta map** in inputs and outputs for tracking

- **Always emit versions.yml** for reproducibility

- **Always include stub block** for fast testing

- **Use resource labels** rather than hardcoding resources

- **Configure via ext.args** not hardcoded parameters

- **One primary tool per module** - compose with subworkflows

- **Use specific container versions** not 'latest'

- **Document container choices** with comments

- **Test with stub mode** before full runs

- **Use consistent naming** - follow meta.id patterns

## Further Reading

- [PATTERNS.md](PATTERNS.md) - Common workflow patterns using modules
- [CONFIGURATION.md](CONFIGURATION.md) - How to configure modules
- [Nextflow Process Docs](https://nextflow.io/docs/latest/process.html)
