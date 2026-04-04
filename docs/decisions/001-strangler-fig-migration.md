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

## Test coverage summary (as of Wave 7 completion)

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
