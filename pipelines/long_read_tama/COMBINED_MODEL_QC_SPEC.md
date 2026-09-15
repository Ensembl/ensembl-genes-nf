# Long-read TAMA combined-model and annotation-QC specification

## Purpose

Produce one non-redundant transcript-model set from all approved long-read
runs, then evaluate that combined set with protein-homology and future
annotation-QC methods. Model classification must be performed only after the
combined model has been generated and all enabled QC statistics have been
joined to it.

The required reduction order is:

```text
reads within one run
    -> TAMA Collapse per run
    -> optional per-run shard merge
    -> TAMA Merge per accession
    -> TAMA Merge across all accessions
    -> one combined transcript set
    -> combined GFF3 and transcript/protein sequences
    -> Diamond and other annotation QC
    -> model classification
```

The final cross-run TAMA merge is the canonical redundancy-reduction step. No
downstream QC process should independently merge, deduplicate, or select a
different transcript set.

## Scope

This specification covers:

- per-run transcript collapse;
- merging contig shards without changing transcript semantics;
- per-accession redundancy reduction;
- cohort-wide redundancy reduction across all approved runs;
- generation of one combined annotation representation;
- sequence derivation for protein-homology searches;
- Diamond as the first model-support QC;
- an extensible interface for future annotation QC;
- classification from the combined model and QC statistics.

It does not define changes to read-input classification, PacBio CCS
generation, minimap2 presets, TAMA parameters, or the scientific thresholds
for accepting/rejecting models. Those are separate decisions and must remain
configuration parameters or documented follow-up work.

## Terminology

- **Run**: one approved ENA accession and its associated long-read input.
- **Run model set**: the result after all reads for one run have been collapsed
  and redundant models from that run have been merged.
- **Combined model set**: the single result after all run model sets in the
  cohort have been merged by TAMA.
- **Model ID**: the identifier assigned to a transcript in the combined model
  set. This is the primary key for QC and classification.
- **Evidence ID**: the originating run, read, or TAMA report identifier used
  for provenance. It is not a substitute for the final Model ID.

## Required workflow topology

### 1. Collapse reads within each run

Each approved run is aligned and passed to `TAMA_COLLAPSE`. The process must
receive all reads for exactly one run and must not mix accessions at this
stage.

The existing default is one TAMA Collapse invocation per run. If
`shard_mode=contig`, each run may first be split by reference contig, but the
shards must be merged back to one run model set before cohort merging.

Required invariant:

```text
one run accession -> one run model set
```

Arbitrary read sharding, overlapping genomic windows, and independent merging
of unrelated chunks are not permitted because they can create duplicate or
truncated models at chunk boundaries.

### 2. Merge contig shards within a run

When contig sharding is enabled:

1. run `TAMA_COLLAPSE` once per contig shard;
2. collect all shard BED files belonging to the same accession;
3. run `TAMA_MERGE` once for that accession;
4. treat that output as the accession’s complete model set.

When sharding is disabled, the output of the single per-run
`TAMA_COLLAPSE` invocation is already the run model set.

The per-run merge must happen before any cross-run merge. It must not be
skipped merely because a run has one shard in a particular execution.

### 3. Merge all runs across the cohort

Collect exactly one run model set per approved run and invoke a dedicated final
`TAMA_MERGE` with all of them in one merge file list.

Required invariant:

```text
all approved runs -> exactly one combined TAMA BED output
```

The combined merge must not be performed independently per tissue, platform,
or accession unless that is explicitly requested as a separate cohort. Any
such partition creates multiple competing canonical model sets.

The final merge file list must be deterministic: sort inputs by run accession
and then by model-set filename before writing it. The sorted list and its
checksum must be published with the results.

## Canonical combined outputs

The final merge stage must publish the following logical products:

| Product | Required content | Consumer |
|---|---|---|
| `combined_models.bed` | Canonical non-redundant TAMA BED12 models | Structural validation, conversion |
| `combined_models.gff3` | One GFF3 representation of the same models | Future annotation QC, reporting |
| `combined_transcripts.fa` | Spliced transcript sequences keyed by Model ID | Sequence QC, ORF prediction |
| `combined_transcripts.faa` | Peptide sequences keyed by Model ID, when available | Diamond blastp, protein QC |
| `combined_model_manifest.tsv` | One row per Model ID with provenance and coordinates | Joins and classification |
| `combined_merge_filelist.tsv` | Exact inputs to the final TAMA merge | Reproducibility |
| `combined_model_stats.tsv` | Counts and structural metrics | QC and release summary |

The BED12 file is the coordinate-level source of truth. The GFF3 and FASTA
files must be generated from that same BED12 file and must preserve a
one-to-one Model ID mapping. They must never be independently regenerated
from source reads or pre-merge models.

## Combined GFF3 contract

The conversion process must emit a valid, deterministic GFF3 containing at
least:

- one `gene` feature per combined gene grouping;
- one `transcript` feature per Model ID;
- one `exon` feature per BED12 block;
- stable `ID` and `Parent` attributes;
- source run/accession provenance;
- chromosome, start, end, strand, and exon rank;
- the original TAMA Model ID in an explicit attribute.

The converter must validate that:

- every BED12 model has a GFF3 transcript;
- every GFF3 transcript maps to one BED12 model;
- exon count and coordinates agree between representations;
- coordinates are valid for the supplied reference FASTA;
- IDs are unique and safe for FASTA and tabular joins.

The initial GFF3 is a transcript/exon annotation. It must not invent CDS
features. CDS/ORF features may be added only by a separately defined,
versioned ORF-prediction or coding-annotation stage.

## Sequence generation

Generate `combined_transcripts.fa` from `combined_models.bed` and the same
reference FASTA used by TAMA. Extraction must be strand-aware and splice BED12
blocks in transcript order.

The sequence-generation process must record:

- reference FASTA checksum;
- reference contig name mapping;
- number of models extracted;
- number of models with invalid coordinates;
- sequence length per Model ID;
- whether the model is complete, partial, or otherwise flagged.

Protein sequence generation is a separate decision because the TAMA model set
does not contain CDS coordinates. The pipeline must support an explicit
protein-generation mode, for example:

1. use supplied CDS/ORF annotations mapped to Model IDs; or
2. run a pinned ORF predictor on `combined_transcripts.fa`; or
3. run a documented six-frame translated search as a screening-only mode.

The selected mode, tool version, parameters, and peptide source must be
recorded in the model manifest. Models without a defensible peptide must not
be silently dropped from the combined transcript set.

## Diamond validation

Diamond is the first consumer of the combined model set. It must run after the
final cross-run merge, not separately on each run’s models for the canonical
classification result.

Inputs:

- `combined_transcripts.faa` or the explicitly selected peptide FASTA;
- a pinned Diamond database built from the selected UniProt/reference protein
  release;
- `combined_model_manifest.tsv`.

The output must be a model-keyed report, for example
`combined_diamond_hits.tsv`, containing at least:

- Model ID;
- query peptide ID;
- best subject ID;
- percent identity;
- aligned length;
- query coverage;
- subject coverage where available;
- e-value;
- bitscore;
- database release/checksum;
- Diamond version and command parameters.

The existing Diamond tabular format can be reused, but the query identifiers
must be guaranteed to map back to final combined Model IDs. A best hit is
evidence, not a classification by itself.

The initial implementation must preserve both hit and no-hit rows by left
joining Diamond results to the complete model manifest. A no-hit model remains
in the combined model set and receives explicit status such as
`NO_PROTEIN_HIT`, `NO_PEPTIDE`, or `NOT_APPLICABLE`.

Diamond thresholds must be configuration parameters and must be benchmarked
against the legacy `HiveBlastRNASeqPep` results before being used for hard
filtering. The comparison must use the same reference database release where
possible and document differences in masking, e-value, target count, and
partial-model handling.

## Future annotation-QC interface

Every future QC method must consume the canonical combined products and emit a
model-keyed result table. It must not create a competing annotation set.

Recommended interface:

```text
QC_METHOD(
    combined_models.bed,
    combined_models.gff3,
    combined_transcripts.fa,
    combined_transcripts.faa when available,
    combined_model_manifest.tsv,
    reference_fasta,
    method-specific databases/configuration
)
    -> <method>_model_stats.tsv
```

Each result table must contain one row per Model ID, including explicit
missing/not-applicable statuses, and must publish tool versions, parameters,
database checksums, and a machine-readable summary.

Potential future consumers include:

- GFF3 structural validation;
- CDS/translation validity;
- reference-annotation concordance;
- BUSCO or lineage-based protein support;
- transcript completeness and length checks;
- splice-junction or alignment support metrics.

## Classification stage

Classification is a downstream join, not part of TAMA Merge or Diamond.

The classifier must left-join:

1. `combined_model_manifest.tsv`;
2. structural model statistics;
3. Diamond statistics;
4. all enabled future-QC statistics;
5. optional provenance/read-support summaries.

It must emit:

- `combined_model_classification.tsv`, one row per final Model ID;
- a classification summary by category;
- a rejected/quarantined-model report, if filtering is enabled;
- the exact rule-set version and threshold configuration.

Classification categories and thresholds must be defined separately from the
workflow wiring. The default policy should be to retain the full combined
model set and attach classifications/flags. Removing models requires an
explicit release policy and a published reason for every removed Model ID.

The classifier must distinguish at least:

- valid structural model;
- protein-supported model;
- protein-unresolved model;
- non-coding or no-ORF model;
- incomplete/low-confidence model;
- failed structural validation.

These categories are not mutually exhaustive until the scientific policy is
approved; the implementation should therefore store independent Boolean or
enumerated evidence fields as well as a derived overall classification.

## Provenance and determinism

The final report bundle must include:

- approved input manifest and checksum;
- one-run model-set manifest;
- final sorted merge file list and checksum;
- TAMA version and parameters for every collapse/merge level;
- reference FASTA checksum;
- combined BED12, GFF3, transcript FASTA, and peptide FASTA checksums;
- Diamond database checksum and release;
- every QC tool version and parameter set;
- classification rule-set version.

Model IDs must remain stable across restart and repeated execution when the
same ordered inputs, reference, TAMA version, and parameters are used.

## Failure and empty-input policy

- A failed run collapse fails the cohort build; it must not silently disappear
  from the final merge.
- A run with zero retained models must produce an explicit empty run model set
  and a reason report.
- The final merge must fail if an approved run is missing from the merge input
  manifest.
- Structural conversion failures fail the combined annotation stage.
- Diamond failure is fatal when Diamond is enabled, but a no-hit result is not
  a pipeline failure.
- Future optional QC methods may report unavailable results only when their
  absence is represented explicitly in the combined QC table.

## Acceptance criteria

The implementation is complete when:

1. every approved run is collapsed independently;
2. every sharded run is merged back to exactly one run model set;
3. exactly one deterministic cross-run TAMA merge produces the canonical
   combined BED12;
4. no downstream process merges or deduplicates models again;
5. the combined GFF3 and transcript FASTA have a validated one-to-one mapping
   to the combined BED12 models;
6. Diamond results can be joined to every final Model ID, including no-hit
   models;
7. future QC methods can be added as independent model-keyed consumers;
8. classification is generated only after all enabled QC results are joined;
9. rerunning with identical inputs and parameters produces identical merge
   inputs, model IDs, and report schemas;
10. a small fixture test proves redundancy is removed both within a run and
    across two runs while distinct isoforms are retained.

## Suggested implementation boundaries

The Nextflow implementation should keep these responsibilities separate:

- `TAMA_COLLAPSE`: reads for one run to run-level BED12;
- `TAMA_MERGE_ACCESSION`: shard/run model sets to one accession model set;
- `TAMA_MERGE_COHORT`: all accession model sets to one combined BED12;
- `BED12_TO_GFF3`: canonical BED12 to combined GFF3;
- `EXTRACT_TRANSCRIPTS`: canonical BED12/reference to transcript FASTA;
- `PREDICT_OR_TRANSLATE_PROTEINS`: explicit peptide-generation strategy;
- `DIAMOND_MODEL_SUPPORT`: combined peptides to model-keyed Diamond stats;
- `ANNOTATION_QC_*`: independent future QC methods;
- `CLASSIFY_COMBINED_MODELS`: final evidence join and classification.

The current `VALIDATE_LONG_READ_MODELS` process should remain as the first
structural check on `combined_models.bed`, but its output should be expanded
eventually from a single count to model-keyed validation results so that
classification can consume it.

