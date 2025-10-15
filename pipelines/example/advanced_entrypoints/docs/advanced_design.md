# Dynamic Entry Point System - Design Document

## Problem Statement

Traditional Nextflow pipelines require developers to manually:
1. Build channels from published files when resuming mid-workflow
2. Track which subworkflows have dependencies on others
3. Validate that required files exist before starting
4. Write repetitive code for each entry point

This creates **high cognitive overhead** and **brittle workflows**.

## Solution: Declarative Entry Point Registry

A centralized registry that declares:
- What each subworkflow **requires** (dependencies)
- What each subworkflow **produces** (outputs)
- How to **automatically resolve** dependencies from published files

## Architecture

### 1. Entry Point Registry (`lib/EntryPoints.groovy`)

```groovy
static final Map SUBWORKFLOWS = [
    'COMBINE_AND_COUNT': [
        description: 'Combine tool outputs and count lines',
        requires: [
            [name: 'tool_a', pattern: '*_A.txt'],
            [name: 'tool_b', pattern: '*_B.txt'],
            [name: 'tool_c', pattern: '*_C.txt'],
            [name: 'tool_d', pattern: '*_D.txt']
        ],
        produces: [
            [name: 'combined', pattern: '*_combined.txt'],
            [name: 'line_counts', pattern: '*_line_count.txt']
        ]
    ]
]
```

**Key insight**: By declaring dependencies explicitly, we can:
- - Auto-detect which entry point to use
- - Validate dependencies before execution
- - Build channels automatically from published files
- - Generate documentation from the registry

### 2. Dependency Validation

```groovy
static Map validateEntryPoint(String baseDir, String entryPoint) {
    def workflow = SUBWORKFLOWS[entryPoint]
    def requirements = workflow.requires

    // Check if required files exist
    def filesExist = checkRequiredFiles(baseDir, requirements)

    if (!filesExist) {
        return [
            valid: false,
            error: "Missing required files: ..."
        ]
    }

    return [valid: true, usePublished: true]
}
```

### 3. Automatic Channel Construction

```groovy
def createChannelFromPublished(String baseDir, List<Map> requirements) {
    def allFiles = []

    requirements.each { req ->
        def pattern = "${baseDir}/${req.name}/${req.pattern}"
        def files = file(pattern)
        allFiles.addAll(files)
    }

    return Channel.fromPath(allFiles)
        .map { file ->
            def sampleId = extractSampleId(file.name)
            [[id: sampleId], file]
        }
}
```

### 4. Main Workflow Orchestration

```groovy
workflow {
    // 1. Determine entry point (auto or manual)
    def entryPoint = params.entry_point == 'auto'
        ? EntryPoints.determineEntryPoint(params.outdir)
        : params.entry_point

    // 2. Validate dependencies
    def validation = EntryPoints.validateEntryPoint(params.outdir, entryPoint)
    if (!validation.valid) {
        error(validation.error)
    }

    // 3. Execute with automatic channel resolution
    switch (entryPoint) {
        case 'COMBINE_AND_COUNT':
            def requirements = SUBWORKFLOWS[entryPoint].requires
            results_ch = createChannelFromPublished(params.outdir, requirements)
            COMBINE_AND_COUNT(results_ch)
            break
        // ... other cases
    }
}
```

## Features

### 1. Automatic Entry Point Detection

```bash
# First run: no outputs exist
nextflow run main.nf --outdir results
# → Auto-detects: FULL workflow

# Second run: tool outputs exist
nextflow run main.nf --outdir results
# → Auto-detects: COMBINE_AND_COUNT (skips tools)
```

### 2. Manual Override

```bash
# Force specific entry point
nextflow run main.nf --entry_point RUN_TOOLS --outdir results
```

### 3. Validation with Clear Errors

```bash
nextflow run main.nf --entry_point COMBINE_AND_COUNT --outdir empty
# → ERROR: Cannot start at entry point 'COMBINE_AND_COUNT':
#    Required published outputs not found in empty
#    Missing directories: tool_a, tool_b, tool_c, tool_d
```

### 4. Self-Documenting

```bash
nextflow run main.nf --help
# → Shows all entry points with their requirements/outputs
```

## Benefits

### For Developers
- **Zero boilerplate**: Declare dependencies once, use everywhere
- **Type safety**: Pattern matching ensures correct files
- **Easy to extend**: Add new entry points to registry
- **Clear errors**: Know exactly what's missing

### For Users
- **Automatic resume**: Pipeline picks up where it left off
- **Explicit control**: Override with `--entry_point` when needed
- **Self-documenting**: `--help` shows all options
- **Fail fast**: Validation before execution

### For Teams
- **Consistent patterns**: All pipelines follow same structure
- **Discoverable**: Entry points documented in code
- **Maintainable**: Changes in one place affect all workflows
- **Testable**: Each entry point can be tested independently

## Extending the System

### Adding a New Entry Point

1. **Add to registry** in `lib/EntryPoints.groovy`:
```groovy
'QUALITY_CONTROL': [
    description: 'Run QC on combined results',
    requires: [
        [name: 'combined', pattern: '*_combined.txt']
    ],
    produces: [
        [name: 'qc_reports', pattern: '*_qc.html']
    ]
]
```

2. **Add case to switch** in `main.nf`:
```groovy
case 'QUALITY_CONTROL':
    def requirements = SUBWORKFLOWS[entryPoint].requires
    combined_ch = createChannelFromPublished(params.outdir, requirements)
    QUALITY_CONTROL(combined_ch)
    break
```

3. **That's it!** Auto-detection, validation, and channel building work automatically.

### Supporting Complex Dependencies

For workflows that need multiple input channels:

```groovy
'ADVANCED_ANALYSIS': [
    description: 'Advanced analysis using multiple inputs',
    requires: [
        [name: 'combined', pattern: '*_combined.txt', channel: 'combined'],
        [name: 'line_counts', pattern: '*_line_count.txt', channel: 'counts']
    ],
    produces: [
        [name: 'analysis', pattern: '*_analysis.json']
    ]
]
```

Then group by channel name when building.

## Future: Nextflow Plugin

This pattern could be generalized into a Nextflow plugin:

### Plugin API Design

```groovy
// In nextflow.config
plugins {
    id 'nf-entrypoints'
}

entrypoints {
    register('COMBINE_AND_COUNT') {
        description 'Combine tool outputs and count lines'

        requires {
            publishDir 'tool_a', pattern: '*_A.txt'
            publishDir 'tool_b', pattern: '*_B.txt'
            publishDir 'tool_c', pattern: '*_C.txt'
            publishDir 'tool_d', pattern: '*_D.txt'
        }

        produces {
            publishDir 'combined', pattern: '*_combined.txt'
            publishDir 'line_counts', pattern: '*_line_count.txt'
        }
    }
}
```

### Plugin Features

- **DSL for declaring entry points**: Clean, readable syntax
- **Automatic validation**: Built into plugin runtime
- **Channel factory**: Generic channel construction
- **CLI generation**: Auto-generate `--help` and parameter validation
- **Tower integration**: Visual entry point selection in UI
- **Schema validation**: Ensure publishDir structure matches requirements
- **Cross-pipeline reuse**: Same plugin works for any pipeline

### Plugin Implementation

```groovy
class EntryPointsPlugin extends BasePlugin {
    void apply(Session session) {
        // Register entrypoint DSL
        session.config.entrypoints.each { name, config ->
            validateEntryPoint(name, config)
            registerEntryPoint(name, config)
        }

        // Add runtime hooks
        session.onWorkflowStart {
            resolveEntryPoint(session)
        }
    }

    void resolveEntryPoint(Session session) {
        def entryPoint = determineEntryPoint(session.config.outdir)
        def validation = validateDependencies(entryPoint)

        if (!validation.valid) {
            throw new EntryPointException(validation.error)
        }

        // Inject channels into workflow
        injectChannels(entryPoint, session)
    }
}
```

## Comparison with Current Approach

### Before (Manual)
```groovy
workflow {
    if (file("${params.outdir}/tool_a").exists() &&
        file("${params.outdir}/tool_b").exists() &&
        file("${params.outdir}/tool_c").exists() &&
        file("${params.outdir}/tool_d").exists()) {

        // Manual channel construction
        tool_a = Channel.fromPath("${params.outdir}/tool_a/*_A.txt")
        tool_b = Channel.fromPath("${params.outdir}/tool_b/*_B.txt")
        // ... more boilerplate

        combined = tool_a.mix(tool_b, tool_c, tool_d)
            .map { file ->
                def id = file.name.replaceAll(/_[A-D]\.txt$/, '')
                [[id: id], file]
            }

        COMBINE_AND_COUNT(combined)
    } else {
        RUN_TOOLS(input_ch)
        COMBINE_AND_COUNT(RUN_TOOLS.out.results)
    }
}
```

### After (Declarative)
```groovy
workflow {
    auto()  // That's it! Everything else is automatic
}
```

## Conclusion

This design demonstrates that **workflow orchestration can be declarative** rather than imperative. By centralizing dependency information in a registry:

1. **Complexity is managed** in one place
2. **Validation is automatic** and comprehensive
3. **Channel construction is generic** and reusable
4. **Documentation is generated** from source of truth

The pattern is **pipeline-agnostic** and could benefit any multi-stage Nextflow workflow. A plugin implementation would make this available to the entire Nextflow ecosystem.
