# ENA Submission Strategy (Ensembl Annotation Evidence)

This document captures decisions and conventions for submitting alignment evidence to ENA via this pipeline.

## Scope
- Submissions are ANALYSIS of type `REFERENCE_ALIGNMENT` pointing to BAM/CRAM aligned to a public assembly.
- One ANALYSIS per annotation / assembly / partial release.
- Each ANALYSIS contains all assessed processed RNA-seq alignment files for that annotation.
- Every submission links to the original data submitter’s source runs with `RUN_REF` and, when known, `SAMPLE_REF`.

## Projects / Studies
- Model: one ENA Project (Study) per assembly + Ensembl partial release (e.g., `prj_GCA_052040795.1-Ensembl-2025-12`).
- The workflow derives and registers this Project automatically before any analyses.
- An umbrella study for Ensembl RNA-seq submissions is expected. The workflow accepts this as `--umbrella_study` and records it as analysis metadata until the official ENA cross-linking model is confirmed.

Rationale: keeps evidence for a given assembly/release discoverable and isolated; avoids alias collisions across releases.

## Idempotency and Retries
- Aliases are stable: re‑submitting the same alias is treated as success (status `EXISTS`).
- Poller detects typical Webin messages like “already exists/alias not unique” and records `EXISTS` with any accession found.
- Safe to re‑run on partial failures; successful accessions are preserved in `accessions.tsv`.

## Linking rules
- RUN_REF: always include in PROD. In TEST, the pipeline omits RUN_REF by default (many runs only exist in PROD).
- SAMPLE_REF: include when the BioSample is known; otherwise omit for cross‑sample evidence.
- EXPERIMENT_REF: optional; use if you prefer experiment‑level links.

## Assembly requirements
- Preferred: `assembly_accession` (GCA/GCF) present in ENA. The XML uses `<ASSEMBLY><STANDARD accession=.../>`.
- If the assembly is not yet in ENA, you can supply `ref_seqs` (comma list) to populate `<SEQUENCE accession=.../>` under `REFERENCE_ALIGNMENT` as an interim reference.

## Manifest shape
The workflow uses a two-step manifest:

- `manifest.tsv`: one row per annotation analysis, with `files_tsv`, `assembly_accession`, `last_geneset_update`, `partial_release_label`, project/study fields, and analysis-level metadata.
- `files.tsv`: one row per BAM/CRAM with local file path, remote filename, run accession, sample accession, and optional header-derived alignment metadata.

Defaults applied by the workflow:
- `analysis_type` coerced to `REFERENCE_ALIGNMENT` (ENA schema requirement).
- In TEST mode, `RUN_REF` omitted unless `omit_run_refs_in_test=false` for the row.

## Operational notes
- Upload via FTP with `lftp`. Aspera is not implemented in the current workflow.
- Webin v2 async queue is used; receipts are polled until success/failure.
- Outputs: XML per analysis, queue responses, and consolidated `accessions.tsv`.

## Open questions
- Final umbrella study accession/alias and whether child project-to-umbrella linkage should be emitted in ENA XML.
- Whether production submissions should use BAMs as-is or CRAMs.
- Whether `SAMD...` sample accessions validate as `SAMPLE_REF`, or need mapping to ENA sample accessions.
- Which registry fields will replace CLI arguments for assembly/annotation metadata.
