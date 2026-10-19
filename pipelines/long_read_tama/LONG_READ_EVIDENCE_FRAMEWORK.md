# Long-read evidence and transcript-confidence framework

## Purpose

The goal of this work is not simply to compare transcript-model counts from
TAMA, StringTie, IsoQuant, FLAIR, Bambu, or other tools. The larger objective
is to establish how long-read transcriptomic data can be used as interpretable
evidence for Ensembl transcript annotation, regardless of whether the reads
come from PacBio or Oxford Nanopore (ONT).

The desired output is a framework that can answer, for an individual
transcript or locus:

1. What parts of the proposed model are directly supported by the reads?
2. How reproducible is that support across libraries, accessions, and
   technologies?
3. Which parts are strongly supported, weakly supported, unresolved, or merely
   not observed?
4. Does a transcript-construction tool make that evidence more interpretable
   than processing the alignments directly?

The framework must distinguish evidence for transcript structure from absence
of evidence. A transcript, intron chain, or UTR must not be rejected merely
because it was not observed in a particular long-read library.

## Two related but distinct prongs

### 1. Replicate and improve the Ensembl long-read implementation

The first prong is engineering and operational. It replicates the useful
functionality of the existing Ensembl long-read pipeline while removing
limitations that make it difficult to use as a general evidence source.

This includes:

- accepting and validating multiple long-read representations, including
  PacBio CCS/HiFi, PacBio processed full-length products, PacBio subreads after
  molecule-aware CCS generation, and ONT FASTQ;
- choosing alignment and preprocessing settings from the observed molecular
  representation rather than trusting filenames or legacy platform metadata;
- using a common reference and alignment contract while retaining technology-
  specific settings where they are necessary;
- preserving native products from multiple transcript-construction tools;
- making large or fragmented assemblies tractable through deterministic
  workload planning and grouped sharding;
- retaining read, junction, end, model, parameter, version, and resource
  provenance;
- making failures, skipped accessions, empty outputs, and partial runs visible;
  and
- producing stable normalized products that can be consumed by downstream
  Ensembl QC and gene-building processes.

This prong establishes that the pipeline is reliable and reusable. It does not
by itself establish the biological value of any backend or data type.

### 2. Use long reads for annotation QC and transcript confidence

The second prong is scientific. It asks what long reads add to annotation
quality control beyond ordinary short-read or alignment summaries.

The central object is not a tool-specific transcript count. It is an evidence
record attached to a candidate transcript and its component features. The
record should support separate confidence assessments for:

- transcript existence;
- intron chain and splice-junction structure;
- transcript start and end coordinates;
- 5' and 3' UTR structure;
- single-exon structure;
- alternative internal exons;
- read and accession reproducibility; and
- conflicts between candidate models or evidence sources.

The reference annotation should be withheld during the primary discovery and
evidence-construction experiment. It can then be introduced as a separate
downstream QC layer. This prevents a tool from appearing successful simply
because it was guided toward the annotation that is later used to assess it.

## The key comparison

The important comparison is between evidence representations, not just between
software names.

### Direct alignment-derived evidence

Process the sorted minimap2 alignments ourselves to derive, at minimum:

- read counts and primary/supplementary alignment counts;
- aligned-length and mapping-quality distributions;
- splice-junction observations;
- complete intron-chain observations per read;
- exact and partial chain support;
- 5' and 3' alignment ends;
- end-position clustering and heterogeneity;
- strand and contig support;
- non-canonical or ambiguous junctions;
- read support across accessions and technologies; and
- evidence that can be traced back to individual reads and alignment records.

This is the baseline. It should remain interpretable even if every model-
construction backend fails or is unavailable.

### Tool-derived candidate models

Run transcript-construction tools on the same reference-free alignment inputs,
retaining their native evidence products. These tools may collapse reads,
resolve ends, merge models, or apply internal filters. Their output is useful
only if each model can be connected back to:

- the supporting reads;
- the observed intron chain;
- the observed transcript ends where available;
- the contributing accession and technology;
- the tool and parameter version; and
- the native support or assignment files.

TAMA, StringTie2/3, tmerge, IsoQuant, FLAIR, and Bambu should therefore be
treated as alternative interpretations of the same alignment evidence. A
backend that produces more models is not automatically better, and a backend
that produces fewer models is not automatically more conservative or more
accurate.

### Hybrid evidence layer

The likely end point is a backend-neutral evidence layer that combines direct
alignment observations with tool-derived model proposals. The model proposal
and the evidence supporting it must remain separate objects:

```text
reads -> alignments -> direct evidence observations
                   -> backend-specific candidate models
                   -> normalized structural identities
                   -> feature-level confidence records
                   -> annotation QC / gene-building decisions
```

The evidence layer should not force all tools to agree. Disagreement is useful
when it identifies structures requiring review or additional evidence.

## What should be measured

### Intron-chain confidence

For each candidate intron chain, record:

- number of reads supporting the exact chain;
- number of independent accessions supporting it;
- number of technologies supporting it;
- number of supporting reads with complete versus partial coverage;
- junction-level support for every intron;
- canonical and non-canonical splice status;
- support for alternative chains sharing the same locus;
- whether the chain is recovered by one or several model-construction tools;
- whether the chain is sensitive to alignment or filtering parameters; and
- the fraction of supporting reads that contain ambiguous, supplementary, or
  low-quality alignments.

Exact repeated support for a complete chain should be treated differently from
many reads that each support only separate fragments of the chain. A chain can
be biologically plausible while still having insufficient direct evidence to
call it confidently.

### UTR confidence

UTR evidence must be assessed separately from coding-region or intron-chain
evidence. Record:

- observed 5' and 3' read-end positions;
- clustering of ends near the proposed transcript boundaries;
- the number of independent libraries and accessions supporting each end;
- end heterogeneity and evidence for incomplete molecules;
- poly(A)- or library-specific end signals where available;
- whether the proposed UTR is supported by complete reads or only by partial
  alignments;
- whether alternative UTR boundaries are observed; and
- whether different tools preserve, extend, or truncate the observed ends.

The absence of a read ending at a proposed UTR boundary is not evidence that
the UTR is invalid. It may reflect transcript abundance, tissue specificity,
library preparation, degradation, incomplete sequencing, or end-calling
behaviour. UTR confidence should therefore be graded by positive support and
reproducibility, not by a binary observed/not-observed rule.

## PacBio and ONT as a sanity check

The original Ensembl long-read implementation was developed primarily around
the more abundant PacBio data available at the time. ONT provides a useful
sanity check of whether the underlying evidence and model-building framework
generalize beyond that historical technology.

The useful ONT question is not simply whether the pipeline accepts ONT FASTQ.
It is whether ONT produces, under an explicit nanopore alignment contract:

- a comparable conversion from aligned reads to candidate transcript models;
- comparable distributions of supported intron-chain lengths and exon counts;
- interpretable read-to-model and junction evidence;
- stable results across accessions or libraries where available; and
- complementary structures rather than only technology-specific artefacts.

Model yield or the ratio of candidate models to aligned reads is a useful
diagnostic, but not a quality metric by itself. It should be reported together
with read support, exact-chain support, end support, multi-accession
reproducibility, model redundancy, and failure/empty-output status.

If ONT and PacBio produce broadly comparable evidence distributions after
accounting for their different error and end profiles, that supports treating
ONT as a useful additional data type rather than excluding it from Ensembl
long-read evidence processing. It would still not prove that every ONT model
is equivalent to a PacBio model.

## Annotation-free primary experiment

The primary experiment should use:

1. a fixed genome reference and documented minimap2 parameters;
2. approved long-read libraries with explicit technology and molecular
   representation;
3. the same alignments for every backend;
4. no reference annotation during candidate-model construction;
5. direct alignment-derived evidence produced independently of the backends;
6. native backend outputs retained alongside normalized structural models;
7. feature-level support records for intron chains and transcript ends; and
8. an explicit completeness and provenance matrix.

The reference annotation should be introduced only after the candidate and
evidence sets are frozen, for secondary analyses such as:

- support for existing Ensembl transcripts;
- evidence for missing or extended UTRs;
- evidence for alternative intron chains;
- detection of annotation models with weak or conflicting long-read support;
- comparison with projection, RNA-seq, homology, and manual review evidence.

This separation makes it possible to distinguish “the tool reconstructed the
annotation” from “the reads independently support the annotation.”

## How to decide whether tools add value

For each backend, compare its models against the direct evidence baseline and
ask:

- Does it recover complete chains that are difficult to identify directly?
- Does it preserve the read and junction evidence needed for interpretation?
- Does it improve end/UTR representation or obscure it?
- Does it reduce redundant models without discarding supported alternatives?
- Does it produce stable results across accessions and technologies?
- Are its failures understandable and recoverable?
- Does it add a useful confidence signal, or only a different model count?

A backend is valuable if it makes the evidence more interpretable, traceable,
or reproducible. It is not valuable merely because it generates more models.

## Expected deliverables

The mature workflow should produce:

- an alignment QC table;
- a read-level splice and end evidence table;
- an intron-chain support table;
- a candidate transcript/model manifest;
- a model-to-read and model-to-junction support table where available;
- a UTR boundary support table;
- a cross-accession and cross-technology reproducibility table;
- a backend comparison table;
- a feature-level confidence classification with reasons; and
- a secondary annotation-QC report linking evidence to Ensembl models.

Confidence should be multidimensional rather than collapsed immediately into a
single score. A downstream Ensembl process may choose to summarize the axes,
but the underlying evidence and uncertainty must remain available.

## Working hypothesis

Long reads are most useful to Ensembl not as a replacement for annotation or
as a single transcript caller, but as a structured source of positive evidence
for transcript architecture. Direct alignment processing provides the
auditable baseline. Model-construction tools may improve discovery and
organization of that evidence, but their contribution must be demonstrated by
better traceability, structural resolution, and reproducibility.

The desired result of this project is therefore a technology-agnostic,
annotation-independent evidence framework that can use PacBio and ONT data,
support confidence in intron chains and UTRs, preserve uncertainty when
features are unobserved, and feed both future gene building and annotation
quality control.
