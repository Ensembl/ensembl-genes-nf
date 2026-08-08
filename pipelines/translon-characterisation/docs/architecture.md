# Architecture and channel contracts

## Design rule

The pipeline's scheduling unit is an attributed instance, represented by one
JSON document and a `meta` map.  `meta` is never mutated; a subworkflow creates
a derived map (`meta + [instance_id: ...]` or `meta + [peptide_id: ...]`) when a
new unit of work is introduced.  Python under `bin/` owns coordinate, frame,
and claim logic. Nextflow only composes files and commodity tools.

## Subworkflows

| Subworkflow | Input contract | Output contract | Why it is a boundary |
|---|---|---|---|
| `INPUT_PREPARATION` | global trusted input files | normalized interval/reference tuples | coordinate convention is fixed once, before all decisions |
| `MSA_REFERENCE_SETUP` | pinned MAF/index URL + SHA-256 manifest | verified chromosome MAF directory | prevents unversioned genome-alignment downloads from entering extraction |
| `SUBSTRATE_TOPOLOGY` | intervals + GENCODE + genome | `[meta(instance), instance.json]`, unhosted calls | emits one transcript-contingent instance per compatible context |
| `TRANSLATION_RECONCILIATION` | instances + TranslonScorer verdicts | reconciled instances | validates the trusted-frame agreement once |
| `INSTANCE_CONTEXT` | reconciled instances | regulatory and expression axis files | transcript-specific, peptide-independent work fans out |
| `CONSTRAINT_AXES` | reconciled instances + alignment/variation references | MSA, PhyloCSF, ORBL, population, GPN axes | preserves independent evolutionary signals and GPU scheduling |
| `PEPTIDE_PREPARATION` | reconciled instances | distinct `(interval, frame)` peptide tuples plus fan-back map | deduplicates peptide work without losing transcript parentage |
| `PEPTIDE_AXES` | distinct peptide tuples | uniqueness, features, structure, detectability axes | peptide work streams once per peptide, not once per isoform |
| `OBSERVED_EVIDENCE` | peptide tuples + controlled observations | instrument-tagged admissible detections | performs uniqueness-aware existence interpretation separately |
| `PRODUCT_CLUSTERING` | collected peptide products/structures | set-level sequence and structural clusters | first intentional global barrier |
| `ADJUDICATION` | collected attributed instances and clusters | typed instance/interval claims and run manifest | second and final global barrier; enforces the coding wall |

## File contracts

All stage records are JSON. A stage writes an `axis` object with exactly one of
`positive`, `null`, `uninformative`, or `unattributable`, a reason, attribution
target, and input/version provenance. Numeric values are supporting fields and
are never used as a cross-axis aggregate.

```
[ meta, path(instance.json) ]
meta = [ id: run_id, interval_id: ..., transcript_id: ..., frame: ... ]

[ meta, path(peptide.json) ]
meta = [ id: run_id, peptide_id: ..., interval_id: ..., frame: ... ]
```

`SUBSTRATE_TOPOLOGY` may emit a soft `translated-no-compatible-transcript`
record. It is carried to adjudication, rather than dropped.

## Module ownership

One module wraps one primary analysis: topology, translation join, regulatory
geometry, expression, MSA extraction, PhyloCSF, ORBL, population constraint,
GPN, peptide derivation, MMseqs2 uniqueness, peptide features, ESMFold,
detectability, observed-evidence join, product clustering, or adjudication.
Each module emits a result channel and `versions.yml`; each has a real and
faithful stub command. Tool arguments, container references, resources, and
optional execution are configured in `nextflow.config`, not in subworkflows.

## Production comparative-genomics set

One alignment is not reused merely to make the axes look consistent. PhyloCSF++
uses the hg38 UCSC multiz100way alignment restricted to the calibrated
58-placental-mammal parameter set. GPN-Star uses its native Cactus alignment:
P243 is the primary recent-constraint axis and M447 is the companion mammalian
axis. V100 is optional. These axes remain separate through adjudication because
they measure constraint over different evolutionary timescales.

The immutable repository revisions, upstream locations, and decision rules are
recorded in `assets/production-reference-contract.json`. Large alignments are
downloaded once into a shared HPC reference cache and verified before compute
jobs start; they are never fetched independently by each scoring task.
