# Human Pangenome Projection

Projects reference gene annotation (e.g. GENCODE on GRCh38) onto target assemblies
(e.g. the ~94 HPRC assemblies) via whole-genome alignment followed by feature-level
coordinate projection, rescue, and validation.

The projection tool (`hpp`) is byte-identical to the original monolithic
implementation it was refactored from (verified on chr20 GRCh38→CHM13); see
`tests/equivalence/`.

## Design

```
MINIMAP2_WGA            swappable aligner front-end -> PAF (cs:Z long tags)
   -> BUILD_SYNTENY     hpp: parse PAF -> blocks, gap fill, sex chromosomes
   -> PROJECT_FEATURES  hpp: coordinate projection, paralog reassign, missed-locus recovery
   -> VALIDATE_MODELS   hpp: initial structural validation
   -> RESCUE_PROJECTIONS hpp: 4 rescue passes + re-validate + protein QC
   -> REFINE_MODELS     hpp: conflict resolution + audit traces
   -> FINALISE          hpp: statistics + GFF3 / stats output
```

The whole-genome aligner is a **separate module** so it can be swapped (minimap2
now; lastz/other for distant species later) without touching the projection
stages — they consume a standard PAF. The hpp stages thread a serialized **state
directory** between processes; running them in sequence reproduces the monolith
exactly.

Each target assembly flows through the DAG independently (per-target parallelism),
which is how the pipeline scales to the full HPRC cohort.

## Inputs

- `--sample_sheet` CSV with header `target_id,target_fasta` (one row per assembly).
- `--ref_fasta` reference genome FASTA (shared).
- `--ref_gff` reference annotation GFF3 (shared).
- `--chromosomes` (optional) comma-separated subset, e.g. `chr20`.

## Run

```bash
cd pipelines/human_pangenome_projection

# Validate the scaffold without running tools
nextflow run main.nf -stub \
    --sample_sheet test/samplesheet.csv \
    --ref_fasta /path/GRCh38.fa \
    --ref_gff   /path/gencode.gff3 \
    --outdir results

# Real run (requires the hpp container — see containers/Dockerfile)
nextflow run main.nf -profile slurm \
    --sample_sheet samples.csv \
    --ref_fasta /path/GRCh38.fa \
    --ref_gff   /path/gencode.gff3 \
    --chromosomes chr20 \
    --outdir /path/results
```

Outputs per target land in `${outdir}/${target_id}/`: `*.mapped.gff3`,
`*.stats.json`, `*.stats.audit.json`, and (if any) `*.stats.removed.tsv`.

## Container

`hpp` is not on Bioconda; build the image from `containers/Dockerfile`, push it to
a registry, and point `--hpp_container` (or `params.hpp_container`) at the URL.
The image bundles minimap2 + samtools because the projection stages shell out to
them for gap-filling and per-gene rescue alignments. The main WGA module uses a
standard minimap2 biocontainer instead.

## The hpp tool

The Python package lives in `hpp/`. It exposes both a monolithic entry point
(`hpp map`) and the six staged commands the Nextflow modules call
(`hpp build-synteny`, `project-features`, `validate-models`, `rescue-projections`,
`refine-models`, `finalise`). Stage state is serialized losslessly via
`hpp/serialize.py` + `hpp/state.py`. Unit + equivalence tests are in `tests/`.
