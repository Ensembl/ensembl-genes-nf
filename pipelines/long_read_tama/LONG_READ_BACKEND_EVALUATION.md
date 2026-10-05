# Long-read candidate-model backend evaluation

## Purpose

This workflow is a scoping benchmark for converting the same long-read
evidence into candidate transcript/gene models for later Ensembl gene
building. It is not a production annotation workflow and a model count is not
by itself a measure of quality.

Every enabled backend must receive the same reference genome, the same
approved long-read libraries, and the same sorted/indexed minimap2 BAMs. Each
backend keeps its own native collapse, discovery, and merge semantics. BED12
is generated only as a common comparison representation after native model
construction; it must never replace a backend's native merge input.

## What is currently implemented

| Backend | Native input/topology | Main value for gene building | Main cost or risk | Scope status |
|---|---|---|---|---|
| TAMA Collapse/Merge | Contig-sharded BAM-derived BED12, then TAMA accession/cohort merge | Explicit long-read collapse and end-aware model reduction; useful baseline | Can be memory-sensitive and sensitive to malformed/empty loci; shard/recovery policy can change evidence | Implemented; native path is the reference baseline |
| StringTie2 | Contig-sharded BAM -> native GTF -> StringTie2 merge at accession/cohort | Fast, familiar GTF output and useful independent assembly baseline | Merge and abundance heuristics are annotation/parameter sensitive; native attributes are not equivalent to evidence support | Implemented; keep native GTF as source of truth |
| StringTie3 | Same topology using the StringTie3 executable/container | Tests the newer StringTie long-read behaviour independently of StringTie2 | Must not be conflated with StringTie2; version and merge semantics need separate accounting | Implemented; compare separately |
| tmerge | Contig-sharded BAM -> read-level exon GTF -> tmerge collapse/merge | Independent read-level splice/exon reconstruction baseline | Extra BAM-to-GTF conversion; not equivalent to feeding collapsed TAMA models; likely sensitive to read-level alignment detail | Implemented; retain read-level and native GTF evidence |
| IsoQuant | Complete sorted/indexed BAM -> annotation-free or reference-guided IsoQuant; no contig-shard merge | Strong candidate for discovery plus read-to-model evidence and assignments; cohort discovery is directly relevant | More memory/storage; data type and mode must be explicit; output filenames and evidence products are version-sensitive | Implemented; annotation-free cohort is primary benchmark mode |
| FLAIR | Complete sorted/indexed BAM per accession -> BAM-derived junctions -> FLAIR transcriptome -> native FLAIR combine | Independent splice-aware transcript discovery with explicit read/isoform mapping; useful complementary candidate set | Python/large intermediate files; combine and filtering parameters can materially change model recall; needs real-resource benchmark | Implemented; annotation-free transcriptome/combine path |
| Bambu | Complete sorted/indexed BAM per accession/cohort -> annotation-free Bambu discovery -> native GTF/RDS | Context-aware de novo discovery with a distinct scoring model and complementary candidate set | R/Bioconductor dependency footprint; quality thresholds and read filtering need calibration; read tracking can be memory-heavy | Implemented as optional annotation-free backend; real target-cohort validation required |

### Relevant tools deliberately kept secondary

TALON is relevant for a later annotation-aware comparison, but it is not an
annotation-free candidate generator in the same sense as the primary matrix:
its normal workflow requires a TALON database initialized from a reference
annotation and requires alignment details such as forward orientation and MD
tags. It should be added as a separate reference-guided/evidence-assignment
experiment once the Ensembl annotation input and alignment contract are fixed.

Bambu is now included in the primary annotation-free matrix because its
official workflow supports de novo discovery from genomic BAMs when no
annotation is supplied. Its reference-guided and quantification modes remain
separate follow-up experiments; do not mix those results with the primary
candidate-model comparison.

The `all` backend mode runs all seven independent candidate-model paths. It does
not make one backend canonical and does not pass one backend's models through
another backend's merge tool.

## Fair comparison contract

The primary benchmark should use:

1. one fixed reference FASTA and index;
2. one fixed set of approved libraries, with PacBio/ONT technology recorded
   from the manifest rather than inferred from filenames;
3. one minimap2 alignment per library, reused by every backend;
4. no existing annotation in the primary discovery comparison;
5. the same cohort definition and accession inclusion/exclusion decisions;
6. native output plus a normalized BED12/GTF view for structural comparison;
7. a complete/partial status and expected-input audit for every backend; and
8. Nextflow trace/timeline/resource data collected under the same executor and
   resource policy.

For a broad exploration run, use `backend_failure_policy=continue` so a
backend-specific failure does not erase the successful methods. Treat the
failed backend as a first-class result: retain its trace/error/status records
and report it as incomplete rather than comparing only the surviving model
counts.

The minimap2 alignment parameters are part of the benchmark contract, not a
backend-specific tuning knob. Record the preset and secondary-alignment
policy (currently `--minimap2_preset splice:hq` and
`--secondary=no`) in every run. If alternative presets, secondary alignments,
or supplementary-alignment handling are explored, rerun the complete backend
matrix against each alignment set; do not compare one backend on one alignment
policy with another backend on a different policy.

Reference-guided modes, where supported, should be a separate secondary
experiment. They answer a different question: assignment/quantification and
correction against an existing annotation, not annotation-free candidate
model discovery.

## Measurements to collect

### Model structure

- number of genes and transcripts, separately for accession and cohort;
- unique intron chains and exact/partial transcript matches between methods;
- single-exon versus multi-exon models;
- transcript length, exon count, intron length, and terminal coordinates;
- strand, chromosome, and boundary distributions;
- fraction of models supported by multiple reads and multiple accessions;
- model overlap and redundancy after coordinate-normalized clustering.

### Evidence and suitability for later gene building

- read-to-model assignments and unassigned reads;
- splice-junction support and non-canonical junction counts;
- 5'/3' end support where the backend provides it;
- reproducibility across accessions and library technologies;
- native evidence needed to distinguish a real isoform from a sequencing or
  alignment artefact;
- whether the output can be consumed without losing transcript IDs,
  attributes, and provenance.

### Computational efficiency

- wall time and CPU time per alignment, accession, shard, and merge stage;
- peak RSS and executor memory requested versus used;
- temporary and published storage footprint;
- parallelism available at accession, contig, and cohort stages;
- restartability after one shard or one accession fails;
- sensitivity to cohort size and transcript/model count.

Counts must always be accompanied by completeness, input-read, and tool
version fields. A backend that silently drops failed shards is not directly
comparable with a complete backend.

## Interpretation for future Ensembl development

No backend should be selected solely because it produces the most models. The
likely future design is a candidate-model/evidence layer that can retain more
than one backend's native products, then apply Ensembl-specific filtering,
projection, homology, and manual/automated review downstream.

The first benchmark should therefore identify:

- a primary candidate generator for stable, evidence-rich models;
- one or more complementary generators that recover structures missed by the
  primary method;
- a common evidence schema for read support, splice support, and provenance;
- a deterministic model-clustering/deduplication stage that is independent of
  any one discovery backend; and
- explicit rules for incomplete or conflicting backend results.

## Remaining implementation needed before biological conclusions

1. Run all backends on a small, fully controlled truth/validation set and on
   the target cohort with real containers; record `trace.tsv`, timeline,
   published native products, and software versions.
2. Extend the backend-neutral comparison report over the canonical BED12
   outputs with structural clusters, read/accession support, and junction/end
   evidence. The pipeline already emits per-backend structural summaries,
   pairwise intron-chain comparison, and a trace-resource summarizer. Do not
   collapse models across backends before reporting each backend separately.
3. Enrich the emitted structural model manifests so each normalized model is
   linked to backend, native model ID, accession, input checksum, tool version,
   read support, junction support, and native support files. The current
   manifests provide the stable structural identity and explicitly mark
   support fields as `UNKNOWN` until backend-specific evidence joins are added.
4. Repeat the real-container smoke on the target executor and confirm exact
   output contracts for pinned IsoQuant and FLAIR versions, especially read
   assignment and combine outputs. IsoQuant requires a writable task-local
   HOME in the current BioContainer; FLAIR requires BAM-derived junctions in
   annotation-free mode.
5. Benchmark parameter sensitivity rather than choosing thresholds from one
   run: TAMA collapse/merge settings, StringTie long-read settings, tmerge
   conversion choices, IsoQuant filters, and FLAIR filtering/combine options.
6. Add biological adjudication against independent evidence or a trusted
   annotation only as a secondary analysis: junction truth, full-length
   controls, RNA-seq/projection/homology support, and known positive/negative
   loci.

Until these measurements exist, the current recommendation is to treat TAMA,
StringTie2/3, tmerge, IsoQuant, FLAIR, and Bambu as complementary
candidate-model generators rather than declare a winner. IsoQuant, FLAIR, and
Bambu provide distinct evidence/scoring views that can inform later model
adjudication; TAMA/StringTie/tmerge provide independent native structure
baselines.
