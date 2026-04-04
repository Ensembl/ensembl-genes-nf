# ADR-001: Strangler-Fig Migration from Perl/eHive to Nextflow

**Date:** 2026-04-03
**Status:** Active

## Context

The Ensembl Genebuild annotation pipeline is implemented in Perl/eHive.
The team is migrating to Nextflow + Prefect incrementally using the
**strangler-fig pattern**: replace one subpipeline at a time while keeping
the existing system running.  The bridge module `HiveRunNextflow` reads an
`output_manifest.json` from each Nextflow pipeline and fans out downstream
eHive jobs.

## Decision

Migrate subpipelines wave by wave.  Each Nextflow pipeline:

1. Lives under `pipelines/<name>/`
2. Follows the **riboseq template conventions** (flat `modules/`, flat
   `subworkflows/`, conditional container, `when` clause, `stub` block,
   `versions.yml` output, slurm/local/conda/docker profiles)
3. Writes `output_manifest.json` for eHive bridge dataflow
4. Has Python unit tests for all `bin/` scripts (pytest, ≥ 80 % coverage)
5. Is **not** connected to an Ensembl core DB — file-based I/O only

## Completed waves

### Wave 1 — Core annotation subpipelines

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `long_read` | minimap2, isoseq3 | long-read FASTQ + genome | BAM, collapsed GFF3 | `f8fcdbe` |
| `repeat_masking` | RepeatMasker, RepeatModeler, RED, TRF, dustmasker | genome FASTA | soft-masked FASTA | `f8fcdbe` |
| `igtr` | genblast (IGTR mode) | IMGT FASTA + genome | ig_gene/tr_gene GFF3 | `94609cc` |
| `genblast_homology` | genblast (7-tier) | protein FASTA batches + genome | genblast_1…7 GFF3 | `58b6cc1` |
| `short_ncrna` | cmsearch (Rfam) + BLASTN (miRBase) | UNMASKED genome | ncRNA GFF3 | `4f07f19` |
| `best_targeted` | exonerate cdna2genome + protein2genome | cDNA/protein FASTA + softmasked genome | best_targeted GFF3 | `ccfcec9` |

### Wave 2 — Assembly and import subpipelines

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `refseq_import` | wget + Python GFF3 parser | GCA/GCF accession | RefSeq GFF3 | `4cc1efa` |
| `load_assembly` | wget + samtools faidx | GCA accession | genome FASTA, synonyms TSV, metadata JSON | `b37ed74` |

### Wave 3 — RNA-seq

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `rnaseq` | STAR + StringTie2 | FASTQ sample sheet + genome | merged rnaseq GFF3 | `c145218` |

### Wave 4 — Projection / LASTZ

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `projection` | LASTZ + axtChain + Python chain projection | softmasked genomes + source GFF3 | projected_transcript GFF3 | `5b8fc76` |

### Wave 5 — Ab initio gene prediction

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `ab_initio` | Augustus | softmasked genome + species model | ab_initio GFF3 | `7f4e7b5` |

### Wave 6 — Geneset consolidation

| Pipeline | Tool(s) | Input | Output | Commit |
|---|---|---|---|---|
| `consolidate` | Python layer-annotation | multiple GFF3 + priority map | consolidated GFF3 | `c745ca2` |

### Wave 7 — eHive PipeConfig glue + local profile fixes

Added eHive `PipeConfig` modules for all 11 completed Nextflow pipelines so
each can be launched by `init_pipeline.pl` without any bespoke Perl wrapper.
Also fixed the `local` profile in the three pre-existing pipelines that had
`singularity { cacheDir = '...' }` instead of `singularity.enabled = false`.

**PipeConfigs added** (`ensembl-analysis` branch `feature/NextflowRunnable`, commits `89933b1b8`, `1ec2316a6`):

| PipeConfig | Nextflow pipeline | Notable params |
|---|---|---|
| `LoadAssembly_conf` | `load_assembly` | `assembly_accession`, `assembly_name` |
| `RnaSeq_conf` | `rnaseq` | `genome_fasta`, `sample_sheet`; optional `star_index`, `annotation_gtf` |
| `BestTargeted_conf` | `best_targeted` | `genome_fasta`; optional `cdna_fasta`, `protein_fasta` |
| `Projection_conf` | `projection` | `source_fasta`, `query_fasta`, `source_gff3`; optional `chain` |
| `AbInitio_conf` | `ab_initio` | `genome_fasta`, `species`; optional `augustus_config_path`, `extrinsic_cfg` |
| `Consolidate_conf` | `consolidate` | `gff3_dir` or `gff3_files` |
| `GenblastHomology_conf` | `genblast_homology` | `genome_fasta`, `uniprot_fasta` |
| `RefseqImport_conf` | `refseq_import` | `assembly_refseq_accession`, `assembly_name`; optional `synonyms_tsv` |
| `ShortNcrna_conf` | `short_ncrna` | `genome_fasta`, `rfam_cm`; optional `mirna_fasta`, `mirna_blast_db` |
| `Igtr_conf` | `igtr` | `genome_fasta`, `igtr_proteins` |
| `LongRead_conf` | `long_read` | `sample_sheet`, `genome_fasta`, `protein_db`; optional `genome_index`; `large_long` RC |
| `RepeatMasking_conf` | `repeat_masking` | `genome_fasta`; optional library paths; `large_long` RC (120 h) |

**Local profile singularity fix** (`ensembl-genes-nf` branch `dev`, commit `30c5aa2`):
`igtr`, `long_read`, `repeat_masking` local profiles now set `singularity.enabled = false` and
`docker.enabled = false` (matching all other pipelines since Wave 2).

### Wave 8 — Orchestration bug-fixes + FullAnnotation_conf

Fixed three critical bugs in `FullAnnotation_conf.pm` and four correctness
issues in the Nextflow pipelines (`ensembl-analysis` commit `f4141a956`,
`ensembl-genes-nf` commit `7febdf5`):

**FullAnnotation_conf restructure** (`ensembl-analysis`):

| Bug | Root cause | Fix |
|---|---|---|
| `RunRepeatMasking` triggered 3× | `LoadAssembly` manifest has 3 outputs; all fanned on ch-2 | `nextflow_dataflow_outputs => 0` for LoadAssembly + RepeatMasking; downstream paths hardcoded from known publishDir convention |
| `RunConsolidate` never triggered | `CollectLayerOutputs` accumulated values but had no funnel/semaphore | `FanAnnotationLayers` now uses eHive `1->A` / `A->1` semaphore; RunConsolidate is the funnel triggered when all 9 layers complete |
| `RunRefseqImport` ran before synonyms_tsv existed | RefseqImport was seeded in parallel with LoadAssembly | Moved to annotation fan (stage 3); synonyms_tsv path hardcoded from LoadAssembly outdir convention |

New architecture: `Seed → LoadAssembly → RepeatMasking → FanAnnotationLayers (1→A→9 pipelines, A→1 RunConsolidate) → ConsumeConsolidatedOutput`

**Nextflow pipeline fixes** (`ensembl-genes-nf`):

| File | Bug | Fix |
|---|---|---|
| `refseq_import/modules/parse_refseq.nf` | `synonyms_tsv ?` always truthy for Nextflow path objects | Check `synonyms_tsv.name != 'NO_FILE'` |
| `repeat_masking/modules/bedtools_maskfasta.nf` | Softmasked FASTA only in work dir; path unpredictable | Added `publishDir "${params.outdir}/genome"` |
| `repeat_masking/bin/write_manifest.py` | Used `p.resolve()` (work-dir path, lost after cleanup) | Now uses `outdir/genome/<name>.softmasked.fa` matching publishDir |
| `consolidate/main.nf` | `gff3_dir/*.gff3` only scanned one level | Changed to `gff3_dir/**/*.gff3` (recursive); repeat GFF3s excluded because MERGE_REPEATS has no publishDir |

**Validation milestone**: `load_assembly` pipeline ran end-to-end locally with real GRCh38.p14 data (NCBI FTP, conda profile on macOS): 4/4 tasks ✔ in 11m 33s, valid `output_manifest.json` written.

## Test coverage summary (as of Wave 8 completion)

| Pipeline | Tests |
|---|---|
| igtr | 40 |
| genblast_homology | 27 |
| short_ncrna | 26 |
| best_targeted | 44 |
| refseq_import | 36 |
| load_assembly | 17 |
| rnaseq | 43 |
| projection | 35 |
| ab_initio | 35 |
| consolidate | 26 |
| long_read | 38 |
| repeat_masking | 12 |
| **Total** | **379** |

## Conventions reference

See `docs/workflows.md`, `docs/MODULES.md`, `docs/PATTERNS.md` for the
module-writing conventions all pipelines must follow.

Key points:
- Flat includes: `include { FOO } from '../modules/foo.nf'`
- Container: `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'singularity URL' : 'docker URL' }"`
- Every process: `when: task.ext.when == null || task.ext.when`
- Every process: `stub:` block + `versions.yml` output
- `withName:` selectors qualified with subworkflow: `'ALIGN:MINIMAP2_ALIGN'`
