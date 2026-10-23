# Technology-aware mixed long-read inputs

## Goal

Allow one species/reference workflow run to use all suitable PacBio and Oxford Nanopore transcript reads, while selecting alignment and model settings from each input's recorded technology and library preparation. Keep the accession, sample, library, archive, and processing history attached to every read set so model outputs and QC remain traceable.

## Current behavior

- `MINIMAP2_ALIGN` already chooses `meta.minimap2_preset` for each manifest row, falling back to `params.minimap2_preset`. Thus PacBio CCS/Kinnex and ordinary ONT cDNA rows can be aligned in one run with different presets.
- The alignment command currently passes only `-ax <preset>` plus a workflow-wide secondary-alignment setting. A preset alone cannot express all required differences, such as minimap2's direct-RNA `-k14 -uf` options.
- IsoQuant and Bambu both receive the workflow-wide `params.isoquant_data_type`. A mixed PacBio/ONT run therefore supplies one type for every BAM, even though the BAM metadata retains the per-row technology.
- `--model_backend` accepts one backend name or `all`; it does not accept a selected list such as `flair,stringtie3,isoquant`. `all` also runs other backends.
- PacBio subread FASTQ is now an explicit alignment route and remains distinguishable from consensus reads. It is aligned as raw observations with CCS marked not applied; Kinnex/MAS-seq segmentation remains an explicit caveat. Subread-BAM-to-CCS candidates use the separate CCS route.

## Manifest and metadata contract

Retain existing classifier, accession, sample, checksums, approval, and audit fields. Add or standardize these row-level fields and validate them against an explicit enum:

| Field | Example values | Purpose |
| --- | --- | --- |
| `technology` | `pacbio`, `ont` | Platform family; drives model input type where required. |
| `read_representation` | `ccs`, `hifi`, `kinnex`, `ont_cdna`, `ont_direct_rna`, `raw_subread` | Distinguishes biological/library/read products that need different preparation or mapping. |
| `preprocessing_state` | `processed`, `needs_ccs`, `ccs_not_applied`, `needs_kinnex_segmentation`, `basecall_required` | Prevents raw/intermediate inputs from being silently treated as final reads. |
| `minimap2_args` | empty or validated option string | Per-row options in addition to a per-row preset. Keep shell-safe construction and reject pipeline-owned output/input flags. |
| `isoquant_data_type` | `pacbio_ccs`, `nanopore` | Explicit IsoQuant type for each eligible input. Do not infer from accession names or file suffixes. |
| `observation_unit` | `consensus_molecule`, `subread_record`, `ont_read` | Makes repeated subread observations auditable and avoids describing them as independent molecules. |

Unknown or conflicting metadata should fail validation or remain quarantined for review; it must not silently fall back to a different technology's parameters. Raw subreads should be preserved as a distinct route. CCS conversion should record tool/container/version and source subread accession. Direct alignment of raw Kinnex subreads is permitted only with the unsegmented state recorded as a technical caveat; transcript-level interpretation should use CCS plus segmentation when those inputs become available.

## Workflow changes

### 1. Per-input minimap2 configuration

Extend alignment metadata and `MINIMAP2_ALIGN` to construct the command from row-level preset and options. Keep a global preset as a backward-compatible fallback for older, approved manifests. Record the resolved preset and arguments in the alignment audit/QC output and provenance.

Initial policy defaults to encode and verify:

- PacBio CCS/HiFi Iso-Seq and processed Kinnex: `splice:hq` (confirm strandedness and Kinnex processing state before use).
- ONT cDNA: `splice`.
- ONT direct RNA: `splice` plus the direct-RNA options recommended by the installed minimap2 version, including `-k14 -uf` where appropriate.
- Raw PacBio subread FASTQ: `splice`, with CCS explicitly marked not applied. Do not label these rows `pacbio_ccs`; keep `observation_unit=subread_record`. CCS conversion remains available only from a Subreads BAM. Kinnex/MAS-seq raw FASTQ is unsegmented and must carry that limitation in its technical notes.

Keep index construction once per species reference and reuse it across the mixed read channel.

### 2. Per-input model typing and cohort grouping

Pass the row-level IsoQuant data type into IsoQuant accession jobs and Bambu accession jobs. For cohort jobs, partition/group BAMs by compatible `isoquant_data_type` and run one cohort job per type; do not combine PacBio and ONT BAMs under one global type. Preserve sample/read-group metadata in each group. If a model backend cannot accept a proposed representation, reject or explicitly skip that row with a recorded reason rather than relabeling it.

For raw subreads, IsoQuant/Bambu eligibility remains a separate policy decision. They must not be passed as consensus by default; the selected manifests leave their IsoQuant type blank and the pipeline reports them as unsupported for those backends.

### 3. Selectable backend sets

Allow an explicit backend list (for example `--model_backends flair,stringtie3,isoquant`) while retaining the existing singular parameter for compatibility. Validate duplicates and unsupported names, then enable only requested subworkflows. Keep `all` as the existing convenience value. Ensure candidate merging/comparison includes exactly the selected outputs and that no unselected backend is run as a side effect.

### 4. Provenance and validation

For every aligned run and model job, retain accession, species/reference, sample/library identifiers, technology, representation, preprocessing state, resolved minimap2 configuration, and model-specific data type. Reports should distinguish missing inputs, unsupported routes, and intentionally skipped accessions. Candidate discovery counts must not be reported as approved production inputs.

## Acceptance criteria

1. One approved species manifest containing PacBio CCS, ONT cDNA, and ONT direct-RNA rows produces per-row alignment commands with the expected settings and one shared reference index.
2. Alignment audit output records the effective settings for every accession, including the fallback case for an older manifest.
3. In a mixed-tech run, IsoQuant and Bambu accession jobs use the correct row-level type; cohort mode creates separate compatible groups and never labels the whole mixed cohort with one type.
4. Selecting `flair,stringtie3,isoquant` runs only those backends; existing single-backend and `all` behavior remains available.
5. Raw-subread FASTQ rows are routed for direct alignment with `ccs_not_applied` metadata; needs-CCS BAMs, needs-Kinnex-segmentation rows, and basecall-required rows remain rejected unless explicitly preprocessed. Their observation/preprocessing metadata survives into outputs.
6. Manifest validation, run inventory, alignment QC, model outputs, and combined candidate manifests preserve a common accession/sample identity for traceability.

## Validation cases to add when implementing

- A mixed PacBio CCS + ONT cDNA + ONT direct-RNA manifest, asserting per-row minimap2 commands and audit records.
- A mixed-tech IsoQuant/Bambu cohort, asserting technology-homogeneous groups and correct `--data_type` values.
- Backend selection for one backend, a selected subset, and `all`, asserting only expected processes are scheduled.
- Negative cases for unknown representation, missing direct-RNA metadata, raw subreads with missing CCS-state metadata, and Kinnex rows with incomplete preprocessing state.
- Backward-compatibility case for an approved manifest using the existing global minimap2 preset and singular backend parameter.

## Current safe operating guidance

Today, one species run can mix supported PacBio consensus and ONT cDNA rows for alignment, with per-row presets. FLAIR or StringTie3 can each be selected as a single backend for that mixed set. Direct-RNA rows need their additional minimap2 options, which the current alignment module does not pass per row. IsoQuant can currently be run on one technology per workflow invocation by setting `--isoquant_data_type` accordingly. Running FLAIR, StringTie3, and IsoQuant together as exactly that subset requires separate invocations today; a mixed PacBio/ONT IsoQuant invocation should wait for per-input typing or be split by technology. Raw-subread FASTQs can be aligned directly under the explicit `raw_subread`/`ccs_not_applied` route; they remain unsupported for IsoQuant and Bambu unless separately validated. Kinnex/MAS-seq segmentation limitations must remain visible in the sample sheet.
