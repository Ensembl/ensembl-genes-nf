# Module Reference

Modules are reusable process definitions. Pipeline-specific modules live under each pipeline, while shared examples and helpers live in the repository-level `modules/` directory.

## Shared Modules

| Module | Purpose |
| --- | --- |
| `modules/minimal_example.nf` | Smallest practical process template. |
| `modules/reference_example.nf` | More complete process reference pattern. |
| `modules/move_to_ftp.nf` | Shared helper for moving published files to an FTP target. |

## Pipeline Modules

| Pipeline | Module directory | Notes |
| --- | --- | --- |
| Example | `pipelines/example/modules/` | Teaching modules for simple and composed workflows. |
| Repeat | `pipelines/repeat/modules/` | Fetching assemblies, managing repeat libraries, and repeat annotation tools. |
| Ribo-seq | `pipelines/riboseq/modules/` | Fetching, QC, trimming, alignment, RiboMetric/RiboWaltz, tracks, and unique-read handling. |
| Translon consensus | `pipelines/translon-consensus/modules/` | BED12 standardisation, consensus reporting, and HTML report generation. |

## Module Checklist

When adding a module:

1. Give outputs stable `emit:` names.
2. Pass sample metadata through as `tuple val(meta), ...` when the module is sample-scoped.
3. Add a process label that maps to the shared resource configuration or define a clear pipeline-specific resource override.
4. Publish only files that users need after the run.
5. Emit `versions.yml` when the underlying tool can report a version.
6. Add or update the pipeline page when the module changes user-visible outputs.
