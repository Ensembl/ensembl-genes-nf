# Ensembl Eukaryotic Annotation Nextflow Pipelines

**Template Repository** - Best practices and patterns for building Nextflow pipelines for eukaryotic genome annotation.

> **Note**: This branch contains templates and documentation. Implementation branches contain the actual production pipelines.

## What's Here

This repository provides templates and examples to help you build Nextflow pipelines:

- **Pipeline templates** - Production-ready patterns and structure
- **Example workflows** - Working examples from simple to complex
- **Module templates** - Reusable process definitions
- **Comprehensive documentation** - Guides, patterns, and best practices
- **Advanced patterns** - Entry point system for complex workflows

## Quick Start

Explore the templates by running the example workflows:

```bash
# Navigate to the example pipeline
cd pipelines/example

# Run the simplest example (single process)
nextflow run workflows/simple_workflow.nf -stub --outdir results

# Run the main workflow (subworkflows chained together)
nextflow run main.nf -stub --outdir results

# Explore the advanced entry point system
cd advanced_entrypoints
nextflow run main.nf --help
```

## Repository Structure

```
ensembl-genes-nf/  (template branch)
├── pipelines/
│   └── example/          # Complete example pipeline showing all patterns
│       ├── workflows/    # Example workflows (start here)
│       ├── subworkflows/ # Subworkflow patterns (parallel, sequential)
│       ├── modules/      # Process templates
│       └── advanced_entrypoints/  # Dynamic entry point system
├── modules/              # Module templates for creating new processes
├── subworkflows/         # Subworkflow templates
├── config/               # Configuration templates
└── docs/                 # Complete documentation and guides
```

## Documentation

**[docs/](docs/)** - Complete documentation

Start with:
- **[docs/QUICK_START.md](docs/QUICK_START.md)** - Step-by-step guide from basics to advanced
- **[docs/PATTERNS.md](docs/PATTERNS.md)** - 12 common Nextflow patterns with examples
- **[pipelines/example/README.md](pipelines/example/README.md)** - Example pipeline overview

## Key Features

### Progressive Examples

The example pipeline demonstrates patterns from simple to complex:

1. **Simple workflows** - Single process execution
2. **Subworkflow composition** - Chaining reusable components
3. **Dynamic entry points** - Automatic workflow resumption

### Best Practices

- Modular design with reusable components
- Consistent metadata (`meta`) handling across processes
- Proper channel operations and data flow
- Version tracking for all tools
- Stub mode for fast testing and development

### Advanced Patterns

- **Entry point registry** - Declarative dependency management
- **Automatic validation** - Verify required inputs before execution
- **Smart resumption** - Detect existing outputs and skip completed steps

## Using This Template

**Learning the patterns?** Start with the [Quick Start Guide](docs/QUICK_START.md) and run the example workflows in `pipelines/example/`.

**Building a new pipeline?**
1. Copy the structure from `pipelines/example/`
2. Use module templates from `modules/` as starting points
3. Follow the patterns in [docs/PATTERNS.md](docs/PATTERNS.md)
4. Reference the example implementations

**Adding to an existing pipeline?** Browse `modules/` and `subworkflows/` for reusable components you can adapt.

## Requirements

- Nextflow ≥ 21.04.0 (DSL2)
- Singularity 

## Private Container Registry

Some implementation pipelines use private images from the EBI Docker registry.
Authenticate Singularity or Apptainer outside Nextflow with an EBI/GitLab
account or token that has registry-read permission:

```bash
export EBI_REGISTRY_USER="${EBI_REGISTRY_USER:-$USER}"
read -r -s EBI_REGISTRY_TOKEN
printf '\n'
printf '%s\n' "$EBI_REGISTRY_TOKEN" | singularity registry login \
  --username "$EBI_REGISTRY_USER" \
  --password-stdin \
  docker://dockerhub.ebi.ac.uk
unset EBI_REGISTRY_TOKEN
```

Check access before launching a pipeline:

```bash
singularity registry list
singularity pull docker://dockerhub.ebi.ac.uk/<project>/<image>:<tag>
```

For CI, use masked `APPTAINER_DOCKER_USERNAME` and
`APPTAINER_DOCKER_PASSWORD` secrets. Never put registry credentials in a
sample sheet, Nextflow parameter, or committed configuration file.

## Branch Structure

- **`template`** (this branch) - Templates, examples, and documentation
- **Implementation branches** - Production pipelines for specific annotation workflows

## Contributing to Templates

Contributions to improve templates and documentation are welcome! Consider:
- Adding new pattern examples
- Improving documentation clarity
- Adding module/subworkflow templates
- Enhancing the entry point system

## Resources

- [Nextflow Documentation](https://nextflow.io/docs/latest/)
- [nf-core Guidelines](https://nf-co.re/docs/guidelines)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)

## License

[Add license information]
