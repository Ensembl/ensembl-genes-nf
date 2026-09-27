# Commodity-tool contracts

The Nextflow layer owns scheduling and a pinned execution environment.  Python
owns conversion to/from the standard typed-axis JSONL contract and rejects
unattributable evidence before adjudication.

| Stage | Intended engine | Container policy | Required result contract |
|---|---|---|---|
| `PEPTIDE_UNIQUENESS` | MMseqs2 | Biocontainers MMseqs2, pin digest | Exact/near-exact identity and coverage, never short-peptide E-value alone |
| `MSA_EXTRACT` / `PHYLOCSF` / `ORBL` | UCSC hg38 multiz100way restricted to 58 placental mammals + PhyloCSF++ | Digest-pinned, freely redistributable images | Splice-aware, frame-specific codons; depth and matched-null provenance |
| `POP_CONSTRAINT` | Python + tabix | Python/htslib image | Expected variants, coverage/AN correction, CDS-overlap attribution |
| `GPN_SCORE` | GPN-Star P243 and M447 | digest-pinned CUDA image, immutable model and alignment revisions | Separate primate and mammal localised start/stop and element deltas |
| `PEPTIDE_FEATURES` | permissively redistributable sequence/domain/disorder tools only | separate digest-pinned images | Length regime and Axis-A attribution for every hit |
| `STRUCTURE` | ESMFold | GPU image, model digest | Five-category output, including `uninterpretable` |
| `PRODUCT_CLUSTERING` | MMseqs2 + Foldseek | separate CPU/GPU images as needed | Sequence clusters, structural clusters, proteome-fragment flags |

No image is silently selected by a floating `latest` tag.  Production profiles
must replace the provisional image references with immutable digests recorded
in `tools_manifest`.

SignalP, TargetP, DeepLoc, TMHMM, and other tools requiring separately accepted
terms or manually downloaded licensed packages are outside the production
contract.  `changeloc` informs the fan-out and result-normalisation design only;
its code and restricted predictor stack are not copied into this pipeline.
