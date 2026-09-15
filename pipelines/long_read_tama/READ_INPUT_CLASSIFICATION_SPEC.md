# Read-input classification and molecule-normalisation specification

## Purpose

Update `pipelines/long_read_tama` so that it selects the processing path from the *observed molecular representation*, rather than from the legacy ENA `instrument_platform` field or a file name in the old manifest. The desired property is one molecular observation per record supplied to minimap2 and TAMA:

- an ONT read is one observation;
- a PacBio CCS/HiFi read is one observation; and
- a PacBio raw-subread ZMW is converted to one CCS read before it becomes an observation.

This work is required because historical source manifests can label a run `PACBIO_SMRT` while the delivered reads are actually ONT, and because a valid FASTQ can contain PacBio raw subreads rather than independent molecules. A gzip-valid FASTQ is therefore not sufficient evidence that it is a valid input for the current minimap2/TAMA path.

This document is an implementation handoff. Implement the milestones in order; do not skip to bulk downloading or production alignment.

## Baseline and problem in the current workflow

| Location | Current behaviour | Required change |
| --- | --- | --- |
| `bin/normalise_manifest.py` | Rejects any manifest row that is not `ENA PACBIO_SMRT`. | Preserve the declared platform as untrusted provenance; never use it as the routing decision. |
| `modules/fastq_dl.nf` | Calls `fastq-dl -a <run>`, requires exactly one `*.fastq.gz`, and verifies it against the historical manifest checksum. | Acquire the selected representation recorded by preflight, with its metadata-derived checksum. It must support FASTQ and PacBio BAM as different inputs. |
| `modules/prepare_fastq.nf` | Performs gzip and summary-stat checks only. | Add structural/header inspection and molecule-granularity checks before a FASTQ is admitted to alignment. |
| `main.nf` | Every downloaded FASTQ goes straight to minimap2 under global `splice:hq`. | Consume only approved classifications and branch by representation before minimap2. |

`FASTQ_DL` currently ignores the URL carried in the normalized manifest and rediscovers FASTQ from the run accession. That was reasonable for a cached FASTQ-only workflow but is not adequate for representation-aware ingestion.

## Non-negotiable decision rules

1. Treat old-manifest platform, ENA platform, SRA platform, submitted-file names, and FASTQ headers as independent evidence. Store all of them.
2. The delivered representation wins over a declared platform for technical classification. A disagreement is still a conflict and must remain visible; do not overwrite the historical value.
3. `CONFLICT`, `UNKNOWN`, `MIXED`, and `PACBIO_SUBREAD_FASTQ_ONLY` are quarantined by default. They must not be passed to minimap2 or TAMA.
4. A raw PacBio subread record is never an annotation observation. Do not select a representative subread, cluster sequences heuristically, or feed raw subreads directly to TAMA.
5. Do not implement a FASTQ-only consensus algorithm in this change. The supported raw-subread route is PacBio subreads BAM -> PacBio CCS -> canonical FASTQ. FASTQ-only raw-subread runs remain explicit review cases until a separately validated policy exists.
6. The production workflow consumes an approved, classification-bearing TSV, not the historical candidate manifest. A reviewer must make any exception explicit in that TSV and the decision must be preserved in outputs.
7. Existing downstream alignment, contig sharding, TAMA Collapse, and TAMA Merge semantics are unchanged after a canonical one-read-per-molecule FASTQ has been produced.

## Evidence to collect per run

Preflight must write machine-readable evidence, even when the run is quarantined. Query ENA's read-run file report using the run accession. The minimum requested fields are:

```text
run_accession
instrument_platform
instrument_model
fastq_ftp
fastq_md5
submitted_ftp
submitted_md5
submitted_format
library_strategy
library_source
library_selection
```

Keep the raw response (or a normalized, lossless copy) in the run's audit directory. Semicolon-separated ENA artifact fields must be expanded into one row per submitted artifact while retaining the run accession and order, so that a BAM, BAI, and PBI cannot be confused with one another.

Where the SRA Toolkit is available, perform a secondary SRA audit and record the reported platform and spot group. This is corroborating evidence, not a hard dependency: an unavailable SRA record/tool must be reported as `unavailable`, not silently interpreted as `ONT` or `PACBIO`. The SRA audit is particularly useful for examples such as `SRR29278220`, where the ENA experiment metadata and the archived read representation disagree.

For a candidate FASTQ, inspect a deterministic sample of the first 1,000 complete records by streaming from the beginning of the selected URL. Do not download a multi-hundred-GB derived FASTQ just to classify it. The probe must:

- require complete, parseable FASTQ records and equal sequence/quality lengths for every sampled record;
- record the first token of each header, the number sampled, malformed-record count, and counts matching each signature below;
- stop only after 1,000 complete records or clean EOF; and
- treat a truncated/failed stream as probe failure, not as an empty sample.

The probe is a classification sample, not a replacement for full-file validation. A fully acquired FASTQ still needs gzip/integrity validation and the same structural validation before it reaches minimap2.

### Header signatures

Match against the first whitespace-delimited header token after `@`.

| Signature | Interpretation | Notes |
| --- | --- | --- |
| `^[^/[:space:]]+/[0-9]+/[0-9]+_[0-9]+$` | PacBio raw subread | The fields are movie / ZMW / polymerase interval. The molecule key is exactly movie + `/` + ZMW. |
| `^[^/[:space:]]+/[0-9]+/ccs$` | PacBio CCS/HiFi | One consensus observation per ZMW, subject to the normal duplicate-ID audit. |
| RFC-4122-like UUID as the complete first token, plus ONT-compatible archive or SRA evidence | probable ONT | UUID alone is not enough to declare ONT. Record ONT fields such as `barcode`, `runid`, and `ch` when present in the rest of a header. |
| More than one high-confidence signature in one sampled file | `MIXED` | Quarantine. |
| No decisive signature | `UNKNOWN` | Quarantine unless an explicit future policy is approved. |

Use an ordinary implementation regex equivalent to the table; do not rely on shell word splitting or an `awk NR % 4` parser for production FASTQ. The illustrative `awk` grouping from `ERR12670540` is valid evidence for a four-line ENA FASTQ, but a production validator must support general valid FASTQ layout and validate it.

### Submitted-artifact signatures

Classify file names case-insensitively after URI decoding. Record the exact matching basename and URI as evidence rather than only a Boolean.

| Submitted artifact | Representation evidence | Preferred input |
| --- | --- | --- |
| `*.subreads.bam` | PacBio raw subreads | the BAM and required PacBio sidecars, then CCS |
| `*.ccs.bam`, `*.hifi_reads.bam`, or a validated equivalent | PacBio consensus reads | BAM -> canonical FASTQ; do not run CCS again |
| a FASTQ with the CCS signature | PacBio consensus reads | FASTQ directly |
| a FASTQ with an ONT signature plus metadata support | ONT reads | FASTQ directly |

Do not guess from a generic `.bam`, `.fastq`, instrument model, or the word `PacBio` in a free-text description. Such runs are `UNKNOWN` unless the evidence above or an approved rule establishes the representation.

## Classification contract

Separate *what the reads are* from *whether they may be processed*. The inspector writes both an audit report and a proposed input manifest.

### `run_classification.tsv`

One row per run; values may be `unknown` or `unavailable`, never blank.

```text
run_accession
tissue
description
manifest_declared_platform
ena_platform
ena_instrument_model
sra_platform
sra_spot_group
selected_artifact_uri
selected_artifact_md5
selected_artifact_basename
submitted_representation
header_representation
header_records_sampled
header_subread_count
header_ccs_count
header_ont_uuid_count
header_malformed_count
classification
confidence
status
reason_codes
proposed_action
review_decision
reviewer
reviewed_at
```

Use these controlled values:

| Field | Allowed values |
| --- | --- |
| `submitted_representation`, `header_representation` | `PACBIO_SUBREAD`, `PACBIO_CCS`, `ONT`, `UNKNOWN`, `MIXED`, `NOT_APPLICABLE` |
| `classification` | `PACBIO_SUBREAD_BAM`, `PACBIO_SUBREAD_FASTQ_ONLY`, `PACBIO_RAW_SUBREAD_FASTQ_UNGROUPABLE`, `PACBIO_CCS_BAM`, `PACBIO_CCS_FASTQ`, `PACBIO_PROCESSED_FASTQ`, `ONT_FASTQ`, `CONFLICT`, `UNKNOWN`, `MIXED` |
| `confidence` | `CONFIRMED`, `PROBABLE`, `INSUFFICIENT` |
| `status` | `READY_FOR_REVIEW`, `QUARANTINED`, `APPROVED`, `REJECTED` |
| `proposed_action` | `RUN_CCS_THEN_ALIGN`, `CANONICALISE_CCS_THEN_ALIGN`, `ALIGN_PACBIO_CCS`, `ALIGN_ONT`, `QUARANTINE` |
| `review_decision` | `PENDING`, `APPROVE`, `REJECT` |

`reason_codes` is a semicolon-delimited, controlled vocabulary, for example `SUBMITTED_SUBREADS_BAM`, `HEADER_SUBREAD`, `HEADER_ONT_UUID`, `ENA_SRA_PLATFORM_MISMATCH`, `MANIFEST_HEADER_MISMATCH`, `FASTQ_PROBE_FAILED`, or `NO_SUPPORTED_ARTIFACT`.

### `run_artifacts.tsv`

Write a separate row per file to avoid losing the relationship between semicolon-expanded URI/MD5 lists:

```text
run_accession  source  uri  basename  md5  format  artifact_role  selected
```

`artifact_role` includes `SUBREADS_BAM`, `CCS_BAM`, `PBI`, `BAI`, `FASTQ`, or `UNKNOWN`. `selected` may only be true for the exact payload selected for the proposed action. The report should also retain unselected artifacts.

### Approved production manifest

Create `approved_run_manifest.tsv` from classification rows with `status=APPROVED` and a compatible `review_decision=APPROVE`. It is the only manifest accepted by production alignment. At minimum it contains:

```text
run_accession
tissue
description
classification
proposed_action
selected_artifact_uri
selected_artifact_md5
selected_artifact_basename
minimap2_preset
expected_header_representation
classification_report_sha256
reviewer
reviewed_at
```

The validator must reject a row that is unapproved, has a missing checksum, has `CONFLICT`/`UNKNOWN`/`MIXED` classification, has an incompatible action, or contains an unsafe accession, URI, checksum, or filename. The reference to the hash makes the approval traceable to a specific audit report.

## Classification and routing matrix

| Evidence outcome | Classification | Default production action | minimap2 preset | Notes |
| --- | --- | --- | --- | --- |
| Submitted `.subreads.bam`; PacBio evidence is consistent | `PACBIO_SUBREAD_BAM` | Download BAM/required sidecars -> CCS -> canonical FASTQ -> align | `splice:hq` | The sole supported raw-subread path. |
| FASTQ headers are movie/ZMW/interval but no usable submitted subreads BAM | `PACBIO_SUBREAD_FASTQ_ONLY` | Quarantine | n/a | Do not choose a pass or manufacture a consensus. |
| Submitted CCS/HiFi BAM; metadata is consistent | `PACBIO_CCS_BAM` | Download BAM -> canonical FASTQ -> align | `splice:hq` | Exclude secondary/supplementary alignments when converting; audit read IDs. |
| FASTQ headers are movie/ZMW/ccs; metadata is consistent | `PACBIO_CCS_FASTQ` | Acquire FASTQ -> validate -> align | `splice:hq` | One consensus record should map to one molecule. |
| ONT-compatible metadata plus sampled UUID/ONT headers | `ONT_FASTQ` | Acquire FASTQ -> validate -> align | `splice` | Do not apply ZMW grouping. Do not infer direct-RNA strand options in this work. |
| Any strong signals disagree (for example ENA PacBio, SRA ONT, UUID headers) | `CONFLICT` | Quarantine pending an explicit reviewed override | n/a | `SRR29278220` is the reference case. Record both values. |
| Conflicting header families in a file | `MIXED` | Quarantine | n/a | Split only under a separately approved, provenance-preserving policy. |
| No supported representation or probe failure | `UNKNOWN` | Quarantine | n/a | Never route by declared platform alone. |

`splice` for ONT and `splice:hq` for CCS are defaults to encode the known read-quality distinction, not proof of a fully optimised alignment policy. Keep a per-row `minimap2_preset` in the approved manifest and validate it against a small allow-list. Do not add automatic `-uf`, direct-RNA, or library-strand inference in this change.

## Workflow design

Implement this as two deliberately separated executions. A production run must not discover, classify, and process an unreviewed accession in one step.

```text
historical candidate manifest
          |
          v
normalise candidate rows (retain declared fields as provenance)
          |
          v
ENA resolver + optional SRA audit + artifact expansion
          |
          +--------------------------+
          |                          |
          v                          v
FASTQ header probe              BAM filename/sidecar inspection
          |                          |
          +------------+-------------+
                       v
         run_classification.tsv + run_artifacts.tsv
                       |
                       v
                 human review/approval
                       |
                       v
              approved_run_manifest.tsv
                       |
             +---------+---------+
             |                   |
             v                   v
 PACBIO_SUBREAD_BAM       CCS/HiFi FASTQ/BAM or ONT FASTQ
 acquire -> CCS -> FASTQ      acquire/canonicalise -> FASTQ
             |                   |
             +---------+---------+
                       v
    full validation + read/molecule audit + technology-specific minimap2
                       |
                       v
        existing sorted BAM -> TAMA Collapse -> existing TAMA Merge
```

Suggested boundaries (names may differ, but preserve the contracts):

1. **Candidate-manifest normalisation:** adapt the current whitespace parser only enough to retain its current run/tissue/description extraction and historical URL, MD5, and platform as declared values. It must no longer require `PACBIO_SMRT`. Prefer a versioned tab-separated manifest for new users; retain the legacy parser only as a migration adapter.
2. **Run resolver/inspector:** a lightweight process or Python tool that writes the two audit TSVs above. It must use a persistent metadata cache keyed by accession and captured response, so reruns do not change an audit merely because the remote service changed. Provide an explicit `--refresh_metadata` option.
3. **Reviewer gate:** a validator that reads the approved TSV and refuses unsafe or non-approved rows before any full-size data acquisition starts. It must never amend a row automatically.
4. **Representation-aware acquisition:** replace the FASTQ-only assumption with content-addressed cache paths that include representation and checksum (for example `.../PACBIO_SUBREAD_BAM/<md5>/<basename>`). Locks, completed download checksum checks, and task-local copies remain mandatory. Verify the ENA-selected artifact's checksum rather than a stale historical FASTQ checksum.
5. **PacBio conversion:**
   - For `PACBIO_SUBREAD_BAM`, call a pinned PacBio CCS tool on the original PacBio BAM and required sidecars, then convert the resulting consensus BAM to canonical FASTQ.
   - For `PACBIO_CCS_BAM`, do not call CCS; convert primary consensus records to canonical FASTQ.
   - In both cases, record input/output read counts, distinct movie/ZMW counts when parseable, tool version, CCS arguments, and a checksum for the canonical FASTQ.
   - First validate the exact CCS tool version, necessary `.pbi`/other sidecars, and command options on a small approved PacBio subread BAM. Do not choose CCS thresholds or command flags from file names alone.
6. **Canonical FASTQ validation:** replace the current gzip-only preparation with a validator that confirms complete records, ID uniqueness, expected header representation, and one-per-molecule expectations. For PacBio CCS, report duplicate movie/ZMW IDs as an error. For ONT, report duplicate read IDs as an error. Emit a `molecule_audit.tsv` consumed by neither TAMA nor minimap2 but published with the run's reports.
7. **Alignment branch:** carry `classification`, `proposed_action`, and `minimap2_preset` in the Nextflow `meta` map into `MINIMAP2_TO_SORTED_BAM`. The module must record the effective classification and preset in `alignment_stats.tsv`. No raw subread file can satisfy its input contract.

Keep the existing per-accession TAMA Collapse boundaries. Do not change `shard_mode`, TAMA parameters, secondary-alignment behaviour, or merge policy as part of this data-integrity change.

## Implementation milestones

### Milestone 1 — inventory only

Build and test the candidate normaliser, ENA resolver, SRA audit adapter, artifact expansion, header probe, classifier, and report writers. It must run without acquiring complete FASTQs or BAMs. Publish a cohort-level summary with counts by classification, confidence, and reason code.

**Exit criteria:** every Bos long-read accession has exactly one `run_classification.tsv` row, every reported submitted artifact has a row in `run_artifacts.tsv`, and no raw-platform conflict is hidden or auto-corrected. Review the report before starting Milestone 2.

### Milestone 2 — representative-path validation

Use small, explicitly approved fixtures and no production-scale bulk run. Validate one example for each supported route:

- ONT FASTQ;
- PacBio CCS/HiFi FASTQ or BAM; and
- PacBio subreads BAM -> CCS -> canonical FASTQ.

`ERR12670540` is evidence that a raw-subread class exists, but do not use its 492-GB derived FASTQ as a routine test fixture. Prefer a small approved subreads BAM or a deliberately constructed PacBio-compatible fixture for the mechanical test. A full real-data CCS test requires a planned storage and compute allocation.

**Exit criteria:** each supported route produces a canonical FASTQ whose validator succeeds; the raw-subread route reduces multiple records from an observed movie/ZMW to one output consensus record; the ONT route preserves individual records without ZMW collapsing; and all checksums/provenance reports are published.

### Milestone 3 — controlled production integration

Add the approved-manifest gate and representation-aware acquisition/conversion to `main.nf`. Preserve the existing TAMA DAG downstream of canonical FASTQ. Run structural Nextflow tests, then a small approved end-to-end cohort.

**Exit criteria:** a `CONFLICT` row fails before download/alignment; an approved ONT row uses `splice`; an approved CCS row uses `splice:hq`; an approved raw-subread-BAM row invokes CCS and never exposes raw subreads to minimap2; and no pre-existing TAMA result contract regresses.

### Milestone 4 — cohort rollout

Run inventory across the intended Bos cohort, resolve all quarantine rows by review (approve with recorded evidence or reject), and run only the reviewed approved manifest. Compare molecule counts and TAMA model output against the legacy run as a scientific validation, not merely a workflow-success check.

**Exit criteria:** all production inputs are traceable to a classification report; no rejected/quarantined run appears in alignment/TAMA outputs; and the report explains every input-size or evidence-count change caused by CCS.

## Tests

Add unit tests for the pure parsing/classification code before Nextflow tests. Use static ENA/SRA response fixtures; tests must not depend on live archives.

Required unit cases:

1. `movie/ZMW/start_end` headers classify as `PACBIO_SUBREAD` and extract the exact movie/ZMW key.
2. `movie/ZMW/ccs` headers classify as `PACBIO_CCS`.
3. UUID plus ONT metadata classifies as `ONT_FASTQ`; UUID without supporting evidence is not automatically ONT.
4. A PacBio-declared/ONT-observed `SRR29278220`-shaped case is `CONFLICT` and proposed for quarantine.
5. A submitted `.subreads.bam` with matching PacBio evidence becomes `PACBIO_SUBREAD_BAM`.
6. A subread FASTQ with no submitted BAM becomes `PACBIO_SUBREAD_FASTQ_ONLY` and is quarantined.
7. Mismatched semicolon-separated artifact lists, malformed URLs/checksums, mixed header signatures, invalid FASTQ records, and duplicate output molecule IDs fail explicitly.
8. Approved-manifest validation rejects every status/classification/action combination that is not allowed.

Required workflow test:

```bash
nextflow run pipelines/long_read_tama/main.nf -stub-run -profile stub \
  --manifest pipelines/long_read_tama/test/manifest.txt \
  --reference_fasta pipelines/long_read_tama/test/reference.fa
```

Update that structural test to use an approved-manifest fixture or provide a separate inspect-mode test and an approved-run-mode test. It must exercise both branches without invoking external metadata services in stub mode.

For non-stub validation, use a tiny vetted fixture per Milestone 2, pin the containers/tool versions, and retain command lines and checksums in the test record. Run Nextflow syntax/config/schema validation and the repository's focused long-read checks after each milestone.

## Explicit non-goals

- Do not retrofit a FASTQ sequence-clustering or heuristic subread-collapse algorithm.
- Do not use ENA/SRA metadata to silently relabel a conflicting run.
- Do not split a mixed run automatically.
- Do not redesign TAMA parameters, contig sharding, or transcript-model filtering while changing source integrity.
- Do not claim a raw-subread BAM can be processed until the exact PacBio CCS tool, sidecars, arguments, and output contract have passed Milestone 2.
- Do not bulk-download the existing 492-GB derived FASTQ merely to inspect headers; streamed header sampling and submitted-artifact discovery exist to avoid that cost.

## Deliverables

The implementation is complete only when it provides:

1. Versioned candidate and approved-manifest schemas plus validators.
2. Per-run `run_classification.tsv`, `run_artifacts.tsv`, raw metadata snapshot, header-probe report, and `molecule_audit.tsv` outputs.
3. An explicit reviewer gate before bulk acquisition.
4. A representation-aware cache and acquisition path.
5. Validated ONT, PacBio CCS/HiFi, and PacBio-subreads-BAM routes.
6. Unit fixtures/tests for all decision branches and a structural Nextflow test.
7. Updated `README.md`, `nextflow_schema.json`, and parameter documentation explaining the two-phase invocation and quarantine policy.
