# Ensembl Genes Nextflow Pipelines

This repository collects Nextflow DSL2 pipelines, reusable modules, and shared configuration for Ensembl gene annotation and related analysis work.

The documentation is organised around the questions a new user usually has first: which pipeline should I run, what inputs does it need, where do results go, and how do I maintain it without breaking the shared conventions.

## Pipeline Catalog

| Pipeline | Use it for | Main entry point |
| --- | --- | --- |
| Example | Learning the repository patterns with a small runnable workflow. | `pipelines/example/main.nf` |
| Repeat annotation | Fetching assemblies, generating or reusing repeat libraries, and running repeat annotation tools. | `pipelines/repeat/main.nf` |
| Ribo-seq | Processing ribosome profiling data from acquisition through alignment, QC, offsets, and genome tracks. | `pipelines/riboseq/main.nf` |
| Translon consensus | Building a consensus report from BED12 ORF predictions produced by multiple tools. | `pipelines/translon-consensus/main.nf` |

```{toctree}
:maxdepth: 2
:caption: Start Here

getting-started
repository-structure
configuration
```

```{toctree}
:maxdepth: 2
:caption: Pipelines

pipelines/index
```

```{toctree}
:maxdepth: 2
:caption: Reference

reference/workflow-patterns
reference/modules
reference/development
generated/index
```
