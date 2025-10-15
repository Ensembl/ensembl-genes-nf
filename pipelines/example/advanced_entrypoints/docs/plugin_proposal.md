# nf-entrypoints Plugin Proposal

## Executive Summary

**Proposal**: Create a Nextflow plugin that provides declarative entry point management with automatic dependency resolution.

**Problem**: Every multi-stage pipeline reinvents the same patterns for:
- Detecting which workflow stage to start from
- Validating that required published files exist
- Building channels from published outputs
- Providing clear errors when dependencies are missing

**Solution**: A plugin that makes entry point management **declarative, automatic, and zero-boilerplate**.

## Why This Should Be a Plugin

### 1. Universal Problem

**Every complex pipeline needs this:**
- nf-core pipelines (e.g., rnaseq has preprocessing → alignment → quantification → QC)
- Custom organizational pipelines
- Research pipelines with iterative analysis stages
- Production pipelines with reprocessing workflows

**Current state**: Each pipeline implements its own ad-hoc solution (or doesn't support entry points at all).

### 2. Code Reusability

**Current implementation**: ~200 lines of Groovy in `lib/EntryPoints.groovy`

**As a plugin**: This code is written **once** and works for **every pipeline**.

**Impact**:
- Thousands of pipelines could benefit
- Consistent behavior across ecosystem
- Tested and maintained centrally

### 3. DSL Integration

Plugins can extend Nextflow's DSL, enabling clean syntax like:

```groovy
// In nextflow.config
plugins {
    id 'nf-entrypoints@1.0.0'
}

entrypoints {
    'ALIGNMENT' {
        description 'Align reads to reference genome'
        requires {
            input 'raw_reads'  // Indicates needs raw input
        }
        produces {
            publishDir 'aligned', pattern: '*.bam'
            publishDir 'aligned', pattern: '*.bai'
        }
    }

    'QUANTIFICATION' {
        description 'Quantify gene expression'
        requires {
            publishDir 'aligned', pattern: '*.bam'
            publishDir 'aligned', pattern: '*.bai'
        }
        produces {
            publishDir 'counts', pattern: '*_counts.txt'
        }
    }

    'QC' {
        description 'Run quality control'
        requires {
            publishDir 'counts', pattern: '*_counts.txt'
        }
        produces {
            publishDir 'qc_reports', pattern: '*.html'
        }
    }
}
```

### 4. Runtime Integration

**Plugins have access to Nextflow internals:**

```groovy
class EntryPointsPlugin extends BasePlugin {
    @Override
    void apply(Session session) {
        // Register extension
        session.config.registerExtension('entrypoints', EntryPointsExtension)

        // Hook into workflow lifecycle
        session.onWorkflowStart { event ->
            resolveEntryPoint(event.workflow)
        }

        // Add CLI parameters automatically
        session.params.entry_point = 'auto'
    }

    void resolveEntryPoint(Workflow workflow) {
        def config = workflow.config.entrypoints
        def requested = workflow.params.entry_point

        // Auto-detect or validate
        def entryPoint = requested == 'auto'
            ? autoDetect(workflow.params.outdir, config)
            : validate(requested, workflow.params.outdir, config)

        // Inject channels into workflow context
        if (entryPoint.usePublished) {
            workflow.binding.setVariable('inputs',
                buildChannelsFromPublished(entryPoint))
        }
    }
}
```

### 5. Better Error Messages

**Plugin can integrate with Nextflow's error reporting:**

```
ERROR: Entry point 'QUANTIFICATION' cannot start
├─ Missing required files in /path/to/outdir
├─ Expected: aligned/*.bam
├─ Expected: aligned/*.bai
│
├─ Available entry points:
│  ├─ ALIGNMENT (requires raw input)
│  └─ FULL (requires raw input)
│
└─ Suggestion: Run ALIGNMENT first or use --entry_point FULL
```

### 6. Tower Integration

**Plugins can integrate with Nextflow Tower:**

- Visual entry point selection in Tower UI
- Show dependency graph of entry points
- Display which files are available/missing
- One-click resume from specific stage

### 7. Schema Validation

**Plugin can validate publishDir structure:**

```groovy
entrypoints {
    validate {
        strictMode true  // Fail if publishDir doesn't match schema
        warnOnExtra true  // Warn about unexpected files
    }
}
```

**Runtime checks:**
- Ensure `publishDir` patterns match `produces` declarations
- Warn if files exist that aren't declared
- Validate file naming conventions

### 8. Testing Support

**Plugin can provide testing utilities:**

```groovy
// In tests/test_entry_points.nf
include { testEntryPoint } from 'plugin/nf-entrypoints'

testEntryPoint('QUANTIFICATION') {
    mockPublishedFiles([
        'aligned/sample1.bam',
        'aligned/sample1.bai'
    ])

    expectedOutputs([
        'counts/sample1_counts.txt'
    ])
}
```

## Plugin Architecture

### Core Components

```
nf-entrypoints/
├── src/main/groovy/
│   ├── EntryPointsPlugin.groovy           # Main plugin class
│   ├── EntryPointsExtension.groovy        # DSL extension
│   ├── EntryPointValidator.groovy         # Dependency validation
│   ├── ChannelBuilder.groovy              # Channel construction
│   ├── EntryPointResolver.groovy          # Auto-detection logic
│   └── TowerIntegration.groovy            # Tower support
├── src/resources/
│   └── META-INF/
│       └── MANIFEST.MF                    # Plugin metadata
└── build.gradle                           # Build configuration
```

### API Design

#### 1. Configuration DSL

```groovy
entrypoints {
    // Global settings
    autoDetect true
    strictValidation true
    baseDir params.outdir

    // Define entry points
    register('STAGE_NAME') {
        description 'Human-readable description'

        requires {
            // Option 1: Published files
            publishDir 'dirname', pattern: '*.ext', meta: true

            // Option 2: Raw input
            input 'channel_name'

            // Option 3: Multiple inputs
            publishDir 'dir1', pattern: '*.txt', as: 'input1'
            publishDir 'dir2', pattern: '*.csv', as: 'input2'
        }

        produces {
            publishDir 'output_dir', pattern: '*.result'
        }

        // Optional: Custom validation
        validate { ctx ->
            // Custom logic
            return ctx.fileCount > 0
        }
    }
}
```

#### 2. Workflow Integration

**Automatic mode:**
```groovy
workflow {
    // Plugin automatically:
    // 1. Detects entry point
    // 2. Validates dependencies
    // 3. Injects channels as 'inputs'

    MY_PROCESS(inputs)
}
```

**Explicit mode:**
```groovy
workflow {
    def ep = entrypoint.resolve(params.entry_point)

    if (ep.usePublished) {
        MY_PROCESS(ep.channels.input1, ep.channels.input2)
    } else {
        RAW_INPUT(raw_data)
        MY_PROCESS(RAW_INPUT.out)
    }
}
```

#### 3. Programmatic API

```groovy
// In main.nf
import nextflow.entrypoints.EntryPoint

workflow {
    def ep = EntryPoint.current()

    log.info "Running entry point: ${ep.name}"
    log.info "Using published files: ${ep.usePublished}"

    if (ep.usePublished) {
        log.info "Found ${ep.files.size()} published files"
    }
}
```

### Features

#### 1. Auto-Detection Algorithm

```groovy
def autoDetect(String baseDir, Map config) {
    // Get all entry points sorted by dependency depth
    def candidates = config.entrypoints
        .findAll { it.value.requires.publishDir }
        .sort { -it.value.requires.size() }

    // Find the deepest entry point with satisfied dependencies
    for (ep in candidates) {
        if (validateDependencies(baseDir, ep)) {
            return ep.key
        }
    }

    // Fall back to entry point with no dependencies
    return config.entrypoints.find { it.value.requires.input }?.key
}
```

#### 2. Channel Building

```groovy
def buildChannels(String baseDir, Map requires) {
    return requires.collectEntries { req ->
        def pattern = "${baseDir}/${req.publishDir}/${req.pattern}"
        def files = Channel.fromPath(pattern)

        if (req.meta) {
            files = files.map { file ->
                def meta = extractMetadata(file, req.metaPattern)
                [meta, file]
            }
        }

        [(req.as ?: 'input'): files]
    }
}
```

#### 3. Validation Engine

```groovy
class EntryPointValidator {
    ValidationResult validate(String baseDir, EntryPoint ep) {
        def result = new ValidationResult()

        ep.requires.each { req ->
            def dir = new File("${baseDir}/${req.publishDir}")

            if (!dir.exists()) {
                result.addError("Missing directory: ${req.publishDir}")
                return
            }

            def files = dir.listFiles { f ->
                f.name.matches(globToRegex(req.pattern))
            }

            if (!files || files.length == 0) {
                result.addError("No files matching: ${req.pattern}")
            } else {
                result.addInfo("Found ${files.length} files in ${req.publishDir}")
            }
        }

        return result
    }
}
```

## Benefits

### For Pipeline Developers

- **Write once, use everywhere** - No reimplementing entry point logic
- **Zero boilerplate** - Just declare dependencies in config
- **Automatic validation** - Can't start with missing dependencies
- **Type safety** - Pattern matching at runtime
- **Easy testing** - Mock published files in tests

### For Pipeline Users

- **Automatic resume** - Pipeline picks up where it left off
- **Clear errors** - Know exactly what's missing
- **Self-documenting** - `--help` shows entry points
- **Flexible** - Auto mode or manual override
- **Tower integration** - Visual entry point selection

### For the Nextflow Ecosystem

- **Standardization** - Consistent entry point patterns
- **Reusability** - Shared plugin across all pipelines
- **Maintenance** - Centralized bug fixes and improvements
- **Innovation** - New features benefit everyone
- **Documentation** - One place to learn entry points

## Implementation Plan

### Phase 1: Core Plugin (MVP)
- [ ] Plugin infrastructure
- [ ] Basic DSL for entry point declaration
- [ ] Dependency validation
- [ ] Auto-detection algorithm
- [ ] Channel building from published files
- [ ] CLI parameter integration

### Phase 2: Enhanced Features
- [ ] Multiple input channel support
- [ ] Custom validation functions
- [ ] Metadata extraction from filenames
- [ ] Entry point visualization (DOT graph)
- [ ] Detailed error messages with suggestions

### Phase 3: Ecosystem Integration
- [ ] Nextflow Tower integration
- [ ] nf-core template integration
- [ ] Schema validation
- [ ] Testing utilities
- [ ] Documentation and examples

### Phase 4: Advanced Features
- [ ] Conditional entry points
- [ ] Dynamic dependency resolution
- [ ] Partial resume (within a stage)
- [ ] Cost estimation per entry point
- [ ] Entry point profiles (dev vs prod)

## Comparison with Alternatives

### Current Approach (Manual)

**Pros:**
- Full control
- No dependencies

**Cons:**
- Repetitive code
- Error-prone
- Hard to maintain
- No standardization
- Poor error messages

### Proposed Plugin

**Pros:**
- Declarative
- Automatic
- Standardized
- Well-tested
- Ecosystem-wide benefits

**Cons:**
- Plugin dependency (minor - plugins are standard in Nextflow)
- Initial learning curve (offset by better docs)

## Success Metrics

**Adoption:**
- Used in 10+ nf-core pipelines within 6 months
- 100+ GitHub stars
- Mentioned in Nextflow documentation

**Impact:**
- Reduce entry point code by 90%
- Decrease entry point bugs by 80%
- Improve user experience scores

**Community:**
- Active contributors
- Regular feature requests
- Integration requests from other tools

## Prior Art

### Similar Patterns in Other Systems

**Snakemake**: Has `rule all` and target rules, but no automatic detection
**WDL**: Has task dependencies but no entry point concept
**CWL**: Has workflow inputs but requires manual specification

**Nextflow advantage**: Runtime introspection allows automatic detection that other systems can't do.

### Nextflow Plugins with Similar Scope

**nf-validation**: Schema validation for parameters
**nf-aggregate**: Channel operations
**nf-prov**: Provenance tracking

**nf-entrypoints would complement these** by focusing on workflow orchestration.

## Conclusion

This plugin addresses a **universal problem** in complex Nextflow pipelines. The proof-of-concept demonstrates:

1. **Feasibility** - Works in current implementation
2. **Value** - Reduces boilerplate, improves UX
3. **Generality** - Applies to any multi-stage pipeline

**Recommendation**: Develop as an official Nextflow plugin, potentially collaborating with nf-core for widespread adoption.

**Next Steps:**
1. Present to Nextflow core team
2. Propose at nf-core community meeting
3. Create plugin repository
4. Implement MVP
5. Test with real pipelines
6. Submit to Nextflow plugin registry

---

**Would this be valuable to you and your team?** I believe it would transform how we build and use multi-stage Nextflow pipelines.
