# ENA Submission Strategy (Ensembl Annotation Evidence)

This document captures decisions and conventions for submitting alignment evidence to ENA via this pipeline.

## Scope
- Submissions are ANALYSIS of type `REFERENCE_ALIGNMENT` pointing to BAM/CRAM aligned to a public assembly.
- One ANALYSIS per run (decision: per‑run, not tissue‑merged).
- Every submission links to the original data submitter’s `RUN_REF` and, when known, `SAMPLE_REF`.

## Projects / Studies
- Model: one ENA Project (Study) per assembly + Ensembl release (e.g., `prj_GCA_123456789.1_Ensembl_110`).
- The workflow derives and registers this Project automatically before any analyses.
- Umbrella Project(s) for a species or program can be created manually and cross‑linked outside this pipeline.

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

## Manifest shape (per‑run)
Required columns (see examples/manifest.tsv for a full header):
- `file_path`, `file_type` (`bam|cram`), `study` (Project), `analysis_alias`, `title`, `description`,
- `run_accessions` (single run per row), `assembly_accession`, `sample_accession` (optional),
- `analysis_attributes` (freeform `key=value;` list; `attr_*` columns also supported).

Defaults applied by the workflow:
- `analysis_type` coerced to `REFERENCE_ALIGNMENT` (ENA schema requirement).
- In TEST mode, `RUN_REF` omitted unless `omit_run_refs_in_test=false` for the row.

## Operational notes
- Upload via FTP by default; Aspera is supported if configured.
- Webin v2 async queue is used; receipts are polled until success/failure.
- Outputs: XML per analysis, queue responses, and consolidated `accessions.tsv`.

## Open questions
- Do we want automated cross‑linking to umbrella projects (using Analysis Links) from within this pipeline?
- Policy for enriching `ANALYSIS_ATTRIBUTES` (e.g., tissue labels, mapping stats): current keys are freeform and encouraged.

