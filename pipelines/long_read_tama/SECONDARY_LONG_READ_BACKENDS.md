# Secondary long-read approaches

The primary exploration matrix is deliberately limited to methods that can
turn the same genome-aligned long-read BAMs into annotation-free candidate
models. Other tools are relevant, but answer a different question or require
an input contract that is not currently shared by the matrix.

## TALON — reference-guided evidence assignment

TALON should be evaluated after the primary matrix, using an explicitly frozen
Ensembl annotation release. Its normal workflow initializes a TALON database
from a reference GTF, then annotates reads against that database. TALON also
requires forward-oriented genomic alignments with MD tags. The current shared
minimap2 alignment does not request `--MD`, so adding TALON to `all` without a
separate alignment contract would make the comparison unfair.

Required implementation:

1. Add a separate `reference_guided` evaluation mode.
2. Add an alignment variant with `--MD` and an explicit forward-orientation
   contract, while retaining the primary BAMs unchanged.
3. Initialize one database per reference/annotation release, never one per
   backend run.
4. Preserve TALON's read annotation table, novelty labels, database, and raw
   GTF as native evidence.
5. Compare TALON's known/novel assignment behaviour to Ensembl projection and
   homology evidence, not directly to annotation-free model counts.

## Reference-guided Bambu

Bambu is already included in the primary matrix in annotation-free mode. A
separate reference-guided Bambu run should be added later with an explicit
Ensembl GTF, `discovery`/`quant` settings, NDR threshold, and read-tracking
policy. Its results must be reported separately from the current de novo
candidate set.

## ESPRESSO, FLAMES, and Mandalorion

These are useful follow-up candidates, but should not be added to the weekend
`all` run without dedicated fixtures and pinned containers:

- ESPRESSO is principally an alignment-correction/transcript-discovery
  workflow with a multi-stage contract; its use of corrected alignments and
  optional short-read support needs to be made explicit.
- FLAMES is especially relevant to ONT and single-cell/full-length workflows;
  its technology and barcode assumptions should not be silently applied to a
  bulk PacBio/ONT matrix.
- Mandalorion is most relevant to full-length PacBio/ONT isoform processing
  and has different read/end-quality assumptions from the shared minimap2 BAM
  contract.

For each, the minimum acceptance test is: same reference and approved reads,
documented alignment requirements, native output retained, normalized BED12
output, read-to-model evidence, and resource measurements. Until those are
available, they belong in a secondary technology-specific experiment rather
than the primary comparison.

## Non-generators

Tools such as Mikado or downstream transcript selectors should not be counted
as independent long-read candidate generators. They are selection/adjudication
layers and should be evaluated after candidate models have been generated.
