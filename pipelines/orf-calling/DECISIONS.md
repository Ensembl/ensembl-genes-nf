# ORF Calling: Decisions & Conventions

Last updated: 2026-03-08

Scope: Pipeline `pipelines/orf-calling` providing uniform wrappers for ORF callers consuming outputs of `pipelines/riboseq`.

Decisions
- Separate pipelines: keep `riboseq` as preprocessing; implement callers in `orf-calling` to consume BAM/offsets/annotation.
- Module structure: each tool implemented as three processes — PREP, RUN, PARSE — one task each; no combined multi-step processes.
- Python over shell: prefer small Python scripts (functional style) in `bin/` for input preparation and output parsing; shell is reserved for trivial file moves only in stubs.
- Standardized output: emit `<sample_id>.orf_calls.tsv` with columns
  `sample_id,tool,chrom,start,end,strand,frame,transcript_id,orf_id,score,pval,qval,extra_json` under `results/orf_calls/<tool>/<sample>/`.
- Inputs: transcriptome BAM is default for eukaryotic tools; some callers can accept genomic BAM — modules will accept either via manifest.
- QC integration: optional pass-lengths from `qc_gate` are looked up by sample id when available (`--qc_pass_lengths_dir`).
- Containers: use minimal base for stubs; real runs will pin Bioconda/BioContainers or custom Dockerfiles per tool. GPU tools (e.g., RiboTIE) will get a dedicated CUDA image and `accelerator` profile.

Style & Conventions
- Nextflow DSL2; `publishDir` per-sample; labels `process_low|medium|high` consistent with `riboseq`.
- No classes in Python; small single-file utilities with pure functions.
- Keep params in `nextflow_schema.json`; validate via `nf-schema`.
- One module file per tool; no cross-tool coupling.

Open items
- Pin exact CLI and parameters for each tool in RUN steps once containers are finalized.
- Add A-site profile writer upstream to reduce BAM aggregation for cohort tools.

