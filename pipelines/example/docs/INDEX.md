# Example Pipeline - Complete Index

This pipeline demonstrates Nextflow best practices from basic to advanced patterns.

## 📚 Documentation Index

### Getting Started
- **[QUICK_START.md](QUICK_START.md)** - Start here! Three levels of examples
- **[workflows/README.md](workflows/README.md)** - Simple workflow examples
- **[workflows/PATTERNS.md](workflows/PATTERNS.md)** - Common Nextflow patterns

### Advanced Topics
- **[README.md](README.md)** - Dynamic entry point system
- **[DESIGN.md](DESIGN.md)** - Entry point system design
- **[PLUGIN_PROPOSAL.md](PLUGIN_PROPOSAL.md)** - Why this should be a plugin

## 🎯 Examples by Complexity

### Level 1: Basic (Start Here)
```
workflows/simple_workflow.nf          One process, minimal example
  └─ Uses: modules/tool_a.nf

workflows/minimal_example.nf          Two subworkflows chained
  ├─ Uses: subworkflows/run_tools.nf
  └─ Uses: subworkflows/combine_and_count.nf
```

### Level 2: Subworkflows
```
workflows/two_subworkflows.nf         Custom subworkflows
  ├─ Uses: subworkflows/minimal_subworkflow_example.nf
  └─ Uses: subworkflows/simple_sequential_example.nf

subworkflows/minimal_subworkflow_example.nf    Parallel pattern
subworkflows/simple_sequential_example.nf      Sequential pattern
```

### Level 3: Advanced
```
main.nf                               Dynamic entry point system
  ├─ Uses: lib/EntryPoints.groovy    Entry point registry
  ├─ Uses: subworkflows/run_tools.nf
  └─ Uses: subworkflows/combine_and_count.nf
```

## 📁 File Organization

```
pipelines/example/
│
├─ 📖 Documentation
│   ├── INDEX.md                 ← You are here
│   ├── QUICK_START.md          ← Start here!
│   ├── README.md               ← Entry point system docs
│   ├── DESIGN.md               ← System design
│   └── PLUGIN_PROPOSAL.md      ← Plugin proposal
│
├─ 🚀 Main Workflows
│   ├── main.nf                 ← Advanced: Dynamic entry points
│   └── workflows/
│       ├── README.md           ← Workflow docs
│       ├── PATTERNS.md         ← Common patterns
│       ├── simple_workflow.nf  ← Level 1: Simplest
│       ├── minimal_example.nf  ← Level 1: Basic chaining
│       └── two_subworkflows.nf ← Level 2: Custom subworkflows
│
├─ 🔄 Subworkflows (Reusable)
│   └── subworkflows/
│       ├── run_tools.nf                      ← Parallel pattern
│       ├── combine_and_count.nf              ← Sequential pattern
│       ├── minimal_subworkflow_example.nf    ← Teaching example
│       └── simple_sequential_example.nf      ← Teaching example
│
├─ ⚙️ Modules (Processes)
│   └── modules/
│       ├── minimal_example.nf   ← Template with TODOs
│       ├── tool_a.nf           ← Example tool
│       ├── tool_b.nf           ← Example tool
│       ├── tool_c.nf           ← Example tool
│       ├── tool_d.nf           ← Example tool
│       ├── combine_outputs.nf  ← Utility process
│       └── count_lines.nf      ← Utility process
│
└─ 🔧 Library (Advanced)
    └── lib/
        └── EntryPoints.groovy  ← Entry point registry

```

## 🎓 Learning Path

### Beginner Track
1. Read [QUICK_START.md](QUICK_START.md)
2. Run `workflows/simple_workflow.nf`
3. Examine `modules/tool_a.nf`
4. Run `workflows/minimal_example.nf`
5. Read [workflows/README.md](workflows/README.md)

### Intermediate Track
6. Run `workflows/two_subworkflows.nf`
7. Examine `subworkflows/minimal_subworkflow_example.nf`
8. Read [workflows/PATTERNS.md](workflows/PATTERNS.md)
9. Try implementing your own subworkflow
10. Experiment with different patterns

### Advanced Track
11. Run `main.nf` with different entry points
12. Read [DESIGN.md](DESIGN.md)
13. Examine `lib/EntryPoints.groovy`
14. Read [PLUGIN_PROPOSAL.md](PLUGIN_PROPOSAL.md)
15. Consider contributing to the plugin!

## 🚀 Quick Commands

```bash
# Level 1: Simple examples
nextflow run workflows/simple_workflow.nf -stub --outdir results
nextflow run workflows/minimal_example.nf -stub --outdir results

# Level 2: Subworkflow examples
nextflow run workflows/two_subworkflows.nf -stub --outdir results

# Level 3: Entry point system
nextflow run main.nf --help
nextflow run main.nf -stub --outdir results
nextflow run main.nf -stub --entry_point RUN_TOOLS --outdir results
```

## 📊 Feature Comparison

| Feature | Simple | Minimal | Two Sub | Main (Advanced) |
|---------|--------|---------|---------|-----------------|
| Processes | 1 | 6 | 4 | 6 |
| Subworkflows | 0 | 2 | 2 | 2 |
| Entry points | No | No | No | **Yes** |
| Auto-detection | No | No | No | **Yes** |
| Validation | No | No | No | **Yes** |
| Lines of code | ~30 | ~25 | ~40 | ~130 |
| Complexity | ⭐ | ⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

## 🎯 Use Cases

### Use `simple_workflow.nf` when:
- Learning Nextflow basics
- Testing a single process
- Quick prototyping

### Use `minimal_example.nf` when:
- Learning subworkflow structure
- Building simple linear pipelines
- Understanding channel passing

### Use `two_subworkflows.nf` when:
- Creating reusable subworkflows
- Understanding parallel vs sequential patterns
- Building modular pipelines

### Use `main.nf` (entry points) when:
- Building complex multi-stage pipelines
- Need automatic resume from failures
- Want explicit entry point control
- Building production pipelines

## 🤝 Contributing

Want to add more examples or improve documentation?

1. Create minimal, focused examples
2. Document clearly with comments
3. Include in this index
4. Test with `-stub` mode
5. Submit PR!

## 💡 Key Concepts Reference

### Workflow Structure
```groovy
workflow {
    input_ch = Channel.of(...)
    SUBWORKFLOW(input_ch)
    NEXT_STEP(SUBWORKFLOW.out.results)
}
```

### Subworkflow Structure
```groovy
workflow MY_SUBWORKFLOW {
    take: input_ch
    main: PROCESS(input_ch)
    emit: results = PROCESS.out.results
}
```

### Process Structure
```groovy
process MY_PROCESS {
    publishDir "${params.outdir}/dirname"
    input: tuple val(meta), path(file)
    output: tuple val(meta), path("output.txt"), emit: results
    script: """
        command --input ${file} --output output.txt
    """
}
```

## 🔗 External Resources

- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core Guidelines](https://nf-co.re/docs/guidelines)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
- [Groovy Documentation](https://groovy-lang.org/documentation.html)

---

**Questions?** Check the docs above or explore the code!
**Ready to start?** → [QUICK_START.md](QUICK_START.md)
