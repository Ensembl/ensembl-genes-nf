# ENA Alignment Metadata Contract

This document defines the metadata contract for `ena_submit`. It follows the
ENA project, sample, and reference-alignment model and is intended to be used
as the implementation and review checklist.

## Object model

- One child `PROJECT` represents one species, assembly, and Ensembl release.
- One `ANALYSIS` represents one submitted BAM or CRAM alignment file.
- Every analysis has a `STUDY_REF`.
- A single-source analysis has its source `SAMPLE_REF`, `RUN_REF`, and, when
  known, `EXPERIMENT_REF`.
- An aggregate analysis has all represented source `RUN_REF` values and all
  unique `SAMPLE_REF` values that can be represented accurately. A synthetic
  BioSample must never be created for a computational aggregate.
- If a large aggregate omits `SAMPLE_REF`, its complete run-to-sample mapping
  must be published as provenance and linked from the analysis.

## Required source metadata

Source metadata is resolved by run accession from the ENA Portal API before XML
generation. The resolver records:

```text
run_accession
sample_accession
secondary_sample_accession
experiment_accession
study_accession
secondary_study_accession
scientific_name
tax_id
library_strategy
library_source
instrument_platform
instrument_model
```

Manifest values may be supplied as expected values, but disagreement with ENA
is a validation error unless an explicit, reviewed override is added later.

## Aggregation scopes

```text
single_run
multi_run_single_sample
multi_sample
```

These are the values currently emitted by the resolver. A future cohort
classification must be added only with an explicit mapping rule and tests.

The scope is written to analysis attributes and provenance. Taxonomy must be
consistent across all source runs in an aggregate.

## XML requirements

The Webin V2 document uses:

```xml
<WEBIN>
  <SUBMISSION_SET>...</SUBMISSION_SET>
  <PROJECT_SET>...</PROJECT_SET>
  <!-- or ANALYSIS_SET -->
</WEBIN>
```

For analyses, `EXPERIMENT_REF` precedes `RUN_REF` as required by the analysis
schema. `REFERENCE_ALIGNMENT` must identify the exact assembly and versioned
reference sequences represented by the submitted file.

Project organism metadata belongs under:

```xml
<SUBMISSION_PROJECT>
  <SEQUENCING_PROJECT/>
  <ORGANISM>
    <TAXON_ID>...</TAXON_ID>
    <SCIENTIFIC_NAME>...</SCIENTIFIC_NAME>
  </ORGANISM>
</SUBMISSION_PROJECT>
```

Analysis attributes are supplementary searchable metadata; authoritative
species and taxon information remains on referenced ENA samples.

## Validation gates

Before production submission, validate that:

- all source run accessions resolve;
- sample, experiment, study, species, and taxon relationships are retained;
- single-source analyses do not lose a known BioSample;
- aggregate taxonomy is consistent;
- each submitted file has an MD5 matching the uploaded file;
- each analysis has exactly one BAM or CRAM;
- each analysis references the intended child study;
- project and analysis XML validate against the current ENA schemas;
- Webin TEST receipts identify the expected object type and accession;
- umbrella child-project relationships are verified after update;
- the exact XML, manifest, resolved metadata, provenance, and receipts are
  retained as release artefacts.

## Official ENA references

- [Webin REST V2 programmatic submission](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/programmatic.html)
- [Submitting read alignments](https://ena-docs.readthedocs.io/en/latest/submit/analyses/read-alignments.html)
- [Registering studies](https://ena-docs.readthedocs.io/en/latest/submit/study/programmatic.html)
- [Umbrella studies](https://ena-docs.readthedocs.io/en/latest/submit/study/umbrella.html)
- [ENA project schema](https://github.com/enasequence/webin-xml/blob/master/src/main/resources/uk/ac/ebi/ena/sra/schema/ENA.project.xsd)
- [ENA analysis schema](https://github.com/enasequence/webin-xml/blob/master/src/main/resources/uk/ac/ebi/ena/sra/schema/SRA.analysis.xsd)
