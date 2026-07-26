# Development Guide

## Add or Change a Pipeline

1. Keep the public entry point at `pipelines/<name>/main.nf`.
2. Put process modules in `pipelines/<name>/modules/`.
3. Put larger compositions in `subworkflows/` or `workflows/`.
4. Include the repository root `nextflow.config`.
5. Maintain `nextflow_schema.json` beside the entry point.
6. Add user-facing documentation under `docs/source/pipelines/`.

## Update Parameter Docs

Parameter reference pages are generated from schema files:

```bash
python docs/scripts/render_schema_docs.py
```

Run this after editing any `nextflow_schema.json`.

## Build Documentation Locally

```bash
python -m pip install -r docs/requirements.txt
python docs/scripts/render_schema_docs.py
sphinx-build -W --keep-going -b html docs/source docs/build/html
```

The GitHub Actions documentation workflow runs the same generation step before building the site.

## Keep Docs Useful

Good pipeline documentation should answer:

- What problem does the pipeline solve?
- What files and parameters are required?
- What command should a user try first?
- What directories and report files should they expect?
- Which parts of the code should maintainers edit for common changes?

Prefer short runnable examples over long prose. When a parameter is defined in a schema, avoid duplicating the full reference by hand.
