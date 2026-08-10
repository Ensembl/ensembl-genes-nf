# Translational ORF caller module plan

This is the implementation and validation checklist for
`feature/translon-analysis`. A tool is not considered enabled for production
until its native command completes as a Nextflow process on the mini chr22
fixture and its native output is converted to the shared TSV and BED12
contract.

## Caller inventory

| Tool | Underlying method | Primary alignment | Native output to capture | Current state | Acceptance test |
|---|---|---|---|---|---|
| RiboCode | frame/periodicity statistical test | transcriptome BAM plus derived transcriptome annotation | ORF table plus native reports | native runner, transcriptome-reference preparation, and parser implemented | mini fixture, raw + TSV + BED12; Nextflow exit 0 |
| Ribotricer | periodicity and ORF translation status | transcriptome BAM | `translating_ORFs.tsv` and reports | native runner and native parser implemented | mini fixture, raw table + parser test + Nextflow exit 0 |
| RiboTaper | triplet-periodicity test | genome BAM plus matched RNA-seq BAM | native ORF/BED output | dropped from the unified suite; ORFquant is the successor genome-BAM route | not in scope |
| ORFquant | R/SaTAnn ORF annotation and quantification | genome BAM | R-generated ORF/quantification tables | custom Bioconda-derived image, R compatibility wrapper, and parser implemented | mini fixture, runner + standardiser exit 0; zero-call native result recorded explicitly |
| Rp-Bp | Bayesian periodicity/model selection | FASTQ plus genome annotation and rRNA/adapter resources | predicted ORF BED12 and Bayes-factor outputs | pinned 3.0.1 runner and explicit FASTQ/resource contract implemented; mini fixture lacks required FASTQ/resources | real mini FASTQ/resource fixture, `prepare-rpbp-genome`, `run-all-rpbp-instances`, output discovery, standardisation |
| iRibo | periodicity/scoring caller | genome BAM plus annotation | `translated_orfs.csv`/GFF3 | source-pinned container now builds and native candidate/profile stages complete through Nextflow; sparse mini fixture fails in upstream `GenerateTranslatome.R` with an empty scrambled statistic | adequate real-data smoke, standardised output |
| ORF-RATER | learned ORF classification/scoring | transcriptome BAM plus BED/model assets | scored ORF table | explicit runner added; requires `orfrater_model` and real smoke | source-pinned container and model, real call, standardised output |
| PRICE | probabilistic translation prediction | cohort of genome BAMs plus reference index | `*.orfs.tsv` and companion model/CIT files | nf-core/gedi index + PRICE modules vendored; runner repairs missing BAM MD/NM tags with samtools and tailored native parser; full chr22 run entered model fitting but did not converge, sparse window fails deterministically | adequate real-data cohort/window, native completion, standardized output |
| RibORF | periodicity/frame scoring | FASTQ + genePred + genome | RibORF score table | explicit runner added; BAM-only input is rejected; real smoke still required | source-pinned container, real call, standardised output |
| Ribo-TISH | frame periodicity test and annotation | genome BAM | prediction table/BED | Bioconda-pinned runner and tailored parser implemented | mini fixture, native table + TSV + BED12; Nextflow exit 0 |
| RiboTIE | learned/model-based ORF scoring | transcriptome BAM + reference/model | prediction CSV/GTF | explicit runner added; source/model image and real smoke still required | source-pinned CPU/GPU-compatible container and model, real call, standardised output |

## Work sequence

1. **Shared contract and fixtures**
   - Keep the real mini chr22 BAMs, indexes, FASTA, GTF, and samplesheet immutable.
   - Add any derived tool-specific config or annotation under a separate fixture
     directory, never by modifying the downloaded inputs.
   - Keep the two-process interface: `RUN_<TOOL>` then `STANDARDISE_<TOOL>`.
   - Prefer published RiboSeq BAMs directly. Put any legacy-format conversion
     in `subworkflows/caller_inputs.nf`, keyed by sample, rather than inside the
     main caller graph.

2. **Periodicity callers**
   - Finish Rp-Bp after RiboCode/Ribotricer.
   - Consume read lengths and P-site offsets from the RiboSeq outputs where the
     caller supports them; do not hard-code production offsets.
   - Validate each caller independently before joining outputs.

3. **R/learned/probabilistic callers**
   - Use one pinned container per tool, with model/data assets pinned separately
     from the workflow code.
   - Fail if the executable, model, or expected native output is absent.
   - Do not create an empty or synthetic result to satisfy a Nextflow output.

4. **Standardisation**
   - Parse the native output in a tool-specific adapter where formats differ.
   - Preserve native scores and diagnostic columns in `extra_json`.
   - Emit valid zero-row files when a real run finds no calls, but only after the
     native tool itself has completed successfully.

5. **Validation**
   - Run every tool with Nextflow 26.04.6 and Docker on the mini fixture.
   - Require process exit code 0, non-empty native output directories or an
     explicitly documented zero-call native result, valid standardized headers,
     valid BED12 formatting, and a versions file.
   - Run the consensus boundary with a local memory profile; the existing
     consensus process must not request 60 GB on a 24 GB machine.
   - Record command, image digest/tag, task work directory, output counts, and
     resource use in `.modules-validation-*` logs.

## Distribution notes

Bioconda publishes container images for RiboCode, Ribotricer, Rp-Bp,
ORFquant, and Ribo-TISH. Rp-Bp is a multi-step pipeline rather than a single
ORF-calling executable, so its module must generate the required configuration
and run the profile/prediction stages explicitly. The other published callers
must not be represented as “real” until their source and model assets can be
pinned and exercised on the fixture.
