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

## Test coverage summary (as of Wave 5 completion)

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
| **Total** | **303** |

(long_read and repeat_masking are convention rewrites; bin scripts are Perl
tools wrapped in containers — Python tests not applicable.)

## Conventions reference

See `docs/workflows.md`, `docs/MODULES.md`, `docs/PATTERNS.md` for the
module-writing conventions all pipelines must follow.

Key points:
- Flat includes: `include { FOO } from '../modules/foo.nf'`
- Container: `"${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'singularity URL' : 'docker URL' }"`
- Every process: `when: task.ext.when == null || task.ext.when`
- Every process: `stub:` block + `versions.yml` output
- `withName:` selectors qualified with subworkflow: `'ALIGN:MINIMAP2_ALIGN'`
