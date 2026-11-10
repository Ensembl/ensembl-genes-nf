# Direct alignment and event support

The normalized BED12 files remain the cross-backend transcript-model contract.
With `--evidence_mode enabled`, the long-read pipeline also publishes direct
alignment observations and model/event joins under `evidence/support/`. This
extends the existing evidence graph; it does not filter or rewrite models.

## Source and execution

`ALIGN_LONG_READS` validates the reference and sorted BAM before extracting
features. `SAMTOOLS_VIEW_ALIGNMENT` reads that BAM once into a transient,
compressed SAM stream. `EXTRACT_ALIGNMENT_FEATURES` streams the SAM rows
through the existing feature outputs and a disk-backed support recorder. The
compressed SAM and temporary SQLite database stay in task work directories;
neither is a published evidence product.

The direct-support filter is primary mapped alignment records with MAPQ at
least `evidence_min_mapq` (default 20), equivalent to
`samtools view -F 2308 -q 20` at the default. Duplicate and QC-fail records
remain in this scientific view and are counted separately in the filter flow.
The legacy feature BED/TSV view retains its previous `-F 3844` flag policy.
All counts are alignment-record counts. Distinct read IDs are a separate
quantity, and molecule counts are unavailable without a validated molecule ID.

Introns use zero-based half-open intervals. Raw left and right aligned ends
are zero-based *base* coordinates (`start` and `end - 1`). The raw alignment
orientation is recorded but is not interpreted as transcript orientation.
Model strand determines the 5′ and 3′ assignment; ambiguous model strand is
retained in `model_support` and omitted from strand-specific event tables.

## Published products

The existing evidence package retains `catalog.json`, `relations.tsv`,
`relations.parquet`, the model/alignment artifact graph, and
`verification.json`. The additive `catalog.direct_support` section points to
`support/recording_manifest.json` and `support/SHA256SUMS`. The verifier
checks every support file against its SHA-256 digest and requires one run
directory per alignment manifest.

| Path below `evidence/support/` | Unit and purpose |
|---|---|
| `per_run/<run_id>/alignment_filter_flow.tsv` | Record counts at each filtering/eligibility stage, plus flagged-record counts. |
| `per_run/<run_id>/alignment_noise.tsv` | MAPQ, clipping and available NM-tag summaries. |
| `per_run/<run_id>/chain_observations.tsv.gz` | Exact ordered genomic intron chain, alignment orientation, records and distinct read IDs. Empty chain marks an unspliced alignment. |
| `per_run/<run_id>/junction_observations.tsv.gz` | Junction, chain position, raw support and alignment orientation. |
| `per_run/<run_id>/end_observations.tsv.gz` | Chain-associated exact left/right end pairs and record counts, including unspliced alignments. |
| `summary/model_index.tsv.gz` | One validated BED12 model, original transcript ID where recoverable, structural coordinates and source BED12 checksum. |
| `summary/model_support.tsv.gz` | Model × run exact/partial chain support, minimum junction support, terminal depth, read-set ID and evidence status. |
| `summary/model_junction_support.tsv.gz` | Model × run × junction support and strand-aware splice motif from the declared FASTA. |
| `summary/model_event_support.tsv.gz` | Model × run × 5′/3′ end, exact end support, same-chain depth, endpoint mode and dispersion, and ±5/±10/±20 bp counts. |
| `summary/chain_event_support.tsv.gz` | Deduplicated exact intron-chain event with model, caller, run, sample, library, read-set and technology recurrence; raw depth is the best matching run count, never a sum across callers. |
| `summary/event_caller_run_support.tsv.gz` | Caller × run × exact event; preserves model multiplicity and conditional depth/fraction. |
| `summary/event_support.tsv.gz` | Deduplicated exact 5′/3′ model endpoint event with separate exact-end and ±10 bp raw support, recurrence, and model/caller counts. |
| `summary/event_support_sensitivity.tsv.gz` | Exact-event counts for tolerances 5/10/20 bp and raw-depth thresholds 1/5/10/20/50. |
| `summary/event_cluster_sensitivity.tsv.gz` | Separate single-linkage nearby-end audit; member raw counts are not added. |
| `summary/event_tier_summary.tsv`, `event_best_depth_summary.tsv`, `chain_length_depth_composition.tsv.gz` | Review and depth-composition summaries with explicit denominators. |

Each TSV has a `schema_version` column; direct tables use
`direct-alignment-support-v1`, joined tables use `model-event-support-v1`.
The adjacent metadata JSON files state the coordinate and filter definitions.
The read-level legacy feature files remain under
`per_run/<run_id>/transcript_features/`.

## Interpretation

An exact chain has all model introns in order. A partial chain is a nonempty
contiguous subset of model introns and is never also counted as exact. These
terms describe splice structure, not full-length molecules. Single-exon
support uses interval-overlapping unspliced alignments and is labeled
`single_exon_overlap`; it is not an empty-chain exact match.
The chain-event table uses `partial_only` when only a contiguous subchain is
observed; that tier does not imply full-chain support.

An exact endpoint count requires an aligned end at the model coordinate. The
raw recurrence fields for endpoint events use ±10 bp support, while the
`exact_*` fields use exact-coordinate ends. Chain depth alone does not make
an endpoint supported. The terminal fraction divides records ending within ±10 bp by records in the
same model/run chain stratum. An event's best high-depth fraction is taken
from a model/run row meeting the 20-record reporting threshold. The threshold
and review tiers are sensitivity settings, not a biological truth or accuracy
claim. Caller recurrence shares the underlying reads and is reported
separately from sample and technology recurrence.

`used_to_build` is assigned to run-scope models paired to their run and to
cohort-scope models paired to the contributing runs. These alignments are not
independent validation. `not_observed` means the run was
evaluated without matching evidence; `not_available` means the assay was not
recorded. Full-length and poly(A) status remain unknown/unavailable unless a
future assay-specific source establishes them. Retrospective graph-only
packages have no direct-support section.

## Validation

```bash
pytest -q pipelines/long-read/test
nextflow lint -o concise pipelines/long-read
nextflow config pipelines/long-read/main.nf >/tmp/long-read-event.config
nextflow run pipelines/long-read/main.nf -stub-run -profile test \
  -params-file pipelines/long-read/test/data/stub_all.params.json \
  --evidence_mode enabled --outdir /tmp/long-read-event-stub \
  -work-dir /tmp/long-read-event-work
```

Stub execution checks wiring and schema shape. A real indexed-BAM fixture
checks the arithmetic; a large native Linux/amd64 run is required to establish
resource and container readiness. The four existing HG002 raw-end summaries
can cross-check endpoint counts but cannot validate full chain/junction
coverage because they were restricted to selected caller chains.
