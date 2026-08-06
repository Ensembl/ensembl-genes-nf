# Configuration Strategy

How configuration files are organized and how to use them in your pipelines.

## Configuration Architecture

This repository uses a **layered configuration approach**:

```
Root nextflow.config          ← Base settings (all pipelines)
    ↓ includes
config/resources.config       ← Standard resource labels
config/shared-module.config   ← Module-specific overrides
    ↓ inherited by
pipelines/*/nextflow.config   ← Pipeline-specific settings
```

### Design Principles

1. **Shared defaults at the root** - Common settings inherited by all pipelines
2. **Reusable resource labels** - Standard resource classes defined once, used everywhere
3. **Module configuration separate** - Tool-specific settings isolated in their own file
4. **Pipeline-specific overrides** - Each pipeline can customize as needed

## Configuration Files

### [nextflow.config](../nextflow.config) - Root Configuration

**Purpose**: Base configuration inherited by all pipelines

**What goes here**:
- Default resource allocations (CPU, memory, time)
- Error handling strategy (retry, maxErrors)
- Shell settings (bash flags)
- Execution reporting (timeline, trace, DAG)
- Container/environment defaults (Singularity, Docker)
- Global parameters (outdir, tracedir)

**When to edit**:
- Changing organization-wide defaults
- Updating shared reporting settings
- Modifying container cache locations

### [config/resources.config](../config/resources.config) - Resource Labels

**Purpose**: Standard resource classes for consistent allocation

**What goes here**:
- Process labels with predefined resource allocations
- Standard categories: light, low, medium, high, high_memory, long, etc.

**How to use**:
```groovy
// In your process definition
process MY_ALIGNMENT {
    label 'process_high'  // Automatically gets high resources

    // Process inherits: cpus=8, memory=32GB, time=12h
}
```

**When to edit**:
- Adjusting resource classes for your compute environment
- Adding new resource categories
- Tuning based on actual usage patterns

### [config/shared-module.config](../config/shared-module.config) - Module Settings

**Purpose**: Configure modules without editing their code

**What goes here**:
- Tool-specific command-line arguments (`ext.args`)
- Conditional execution logic (`ext.when`)
- Multiple command configurations (`ext.args2`, `ext.args3`)

**How to use**:
```groovy
// In config file
process {
    withName: 'BLAST' {
        ext.args = '-evalue 1e-10'
        ext.when = { !params.skip_blast }
    }
}

// Module reads these values
process BLAST {
    script:
    def args = task.ext.args ?: ''
    """
    blast ${args} query.fa database
    """
}
```

**When to edit**:
- Configuring tool parameters without touching module code
- Adding conditional execution logic
- Setting different arguments per environment

### Pipeline-Specific Configuration

**Purpose**: Override settings for individual pipelines

**Location**: `pipelines/my-pipeline/nextflow.config`

**What goes here**:
- Pipeline-specific parameters
- Process-specific resource overrides
- Custom profiles for this pipeline
- publishDir configurations

**Example**:
```groovy
// pipelines/annotation/nextflow.config

// Pipeline-specific params
params {
    genome = null
    annotation_db = '/path/to/db'
    min_protein_length = 50
}

// Override resources for specific heavy process
process {
    withName: 'LARGE_GENOME_ALIGNMENT' {
        cpus = 32
        memory = '256.GB'
        time = '72.h'
    }
}
```

## Usage Patterns

### Pattern 1: Use Resource Labels (Preferred)

```groovy
// - Good: Uses standard resource label
process ALIGN {
    label 'process_medium'

    // Automatically gets: cpus=4, memory=16GB, time=8h
}
```

### Pattern 2: Override When Needed

```groovy
// In pipeline's nextflow.config
process {
    withName: 'ALIGN' {
        // Override just what you need
        memory = '64.GB'  // More memory, keeps cpus=4 from label
    }
}
```

### Pattern 3: Configure via ext.args

```groovy
// In config/shared-module.config or pipeline config
process {
    withName: 'TOOL' {
        ext.args = {
            meta.paired_end ?
                '--paired' :
                '--single'
        }
    }
}

// Module stays generic
process TOOL {
    script:
    def args = task.ext.args ?: ''
    """
    tool ${args} input.fq
    """
}
```

## When to Edit What

### Editing Root Config

**Edit [nextflow.config](../nextflow.config) when**:
- Setting organization-wide defaults that all pipelines should inherit
- Configuring shared infrastructure (container cache, work directory)
- Updating execution reporting settings

**Don't edit for**:
- Pipeline-specific parameters
- Process-specific resources
- Tool arguments

### Editing Resource Config

**Edit [config/resources.config](../config/resources.config) when**:
- Tuning resource classes based on your compute environment
- Adding new resource categories for common use cases
- Adjusting allocations based on actual usage data

**Don't edit for**:
- One-off process overrides (use pipeline config instead)

### Editing Module Config

**Edit [config/shared-module.config](../config/shared-module.config) when**:
- Configuring modules used across multiple pipelines
- Setting tool parameters that should be consistent everywhere
- Adding conditional logic that applies broadly

**Don't edit for**:
- Pipeline-specific tool configurations (use pipeline config instead)

### Creating Pipeline Config

**Create `pipelines/my-pipeline/nextflow.config` for**:
- Pipeline-specific parameters
- Process overrides unique to this pipeline
- Pipeline-specific profiles
- Custom publishDir settings

## Best Practices

### Configuration Hierarchy

```
Root defaults
  → Resource labels
    → Module configs
      → Pipeline overrides
        → Command-line params (highest priority)
```

**Key principle**: More specific settings override more general ones

### Do's

- **Use resource labels** - Start with standard labels, override only when necessary

- **Keep modules generic** - Use `ext.args` for configuration, not hardcoded values

- **Document in config files** - Explain why settings exist with inline comments

- **Use meaningful parameter names** - `min_protein_length` not `min_len`

- **Provide sensible defaults** - Pipelines should run with zero configuration

### Don'ts

- **Don't hardcode in modules** - Resources, paths, and arguments belong in config

- **Don't edit root config for one pipeline** - Use pipeline-specific config instead

- **Don't duplicate settings** - If it's in a resource label, don't repeat in process

- **Don't commit local paths** - Use parameters or relative paths

## Configuration Inheritance Example

```groovy
// Root: nextflow.config
process.cpus = 1                    // Default: 1 CPU

// Included: config/resources.config
withLabel:process_medium.cpus = 4   // Label: 4 CPUs

// Pipeline: pipelines/x/nextflow.config
withName:HEAVY_TOOL.cpus = 16       // Override: 16 CPUs

// Result: HEAVY_TOOL gets 16 CPUs (most specific wins)
```

## Further Reading

- [Nextflow Configuration](https://nextflow.io/docs/latest/config.html) - Official docs
- [Process Directives](https://nextflow.io/docs/latest/process.html#directives) - All available settings
