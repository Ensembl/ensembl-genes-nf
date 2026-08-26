# HPC setup for the full translon-analysis suite

This workflow follows the repository convention used by `pipelines/riboseq`:
SLURM for scheduling and Singularity/Apptainer for containers. The supported
Nextflow version is 26.04.6 and Java 17 or 21 is required.

## One-time environment

```bash
module load nextflow
module load apptainer        # or: module load singularity
module load java/21          # Java 17+ is also supported

export NXF_APPTAINER_CACHEDIR=/hps/nobackup/flicek/ensembl/genebuild/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/hps/nobackup/flicek/ensembl/genebuild/singularity_cache
mkdir -p "$NXF_APPTAINER_CACHEDIR"

nextflow -version             # must report 26.04.6
apptainer --version           # or singularity --version
java -version                 # must report 17 or 21
```

Do not build containers in `$HOME` or in the work directory. Pull public
Bioconda images into the shared cache and put custom images in a project
container directory under `/hps/nobackup`.

## Required inputs

The full suite needs:

- the RiboSeq output tree or generated samplesheet;
- genome FASTA and GTF;
- proteome FASTA if characterisation is enabled;
- QC-selected offsets for RibORF;
- ribo-seq FASTQ, rRNA FASTA and adapter FASTA for Rp-Bp;
- an ORF-RATER model directory containing `orfratings.h5`,
  `metagene.txt` and `offsets.txt`;
- a CUDA-capable RiboTIE image if RiboTIE is selected; the current runner
  fine-tunes and predicts from the supplied Ribo-seq BAMs unless additional
  arguments are supplied;
- custom containers for iRibo, ORF-RATER, RibORF, RiboTIE, ORFquant and
  GEDI/PRICE.

The input paths must be visible from compute nodes. Prefer absolute paths under
`/hps/nobackup` and avoid node-local `/tmp` paths for references or models.

## Container preparation

Public images can be pulled directly by Apptainer, for example:

```bash
apptainer pull "$NXF_APPTAINER_CACHEDIR/ribocode.sif" \
  docker://quay.io/biocontainers/ribocode:1.2.15--pyhdc42f0e_1
apptainer pull "$NXF_APPTAINER_CACHEDIR/ribotricer.sif" \
  docker://quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0
apptainer pull "$NXF_APPTAINER_CACHEDIR/ribotish.sif" \
  docker://quay.io/biocontainers/ribotish:0.2.8--pyhdfd78af_0
apptainer pull "$NXF_APPTAINER_CACHEDIR/rpbp.sif" \
  docker://quay.io/biocontainers/rpbp:3.0.1--py310h30d9df9_0
apptainer pull "$NXF_APPTAINER_CACHEDIR/samtools.sif" \
  docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0
```

Build the custom images from the Dockerfiles in `containers/` on a machine
with Docker/BuildKit, push them to an internal registry, and pull them as SIFs
on HPC. Apptainer can also build directly from a Dockerfile where site policy
allows it, but registry-backed immutable tags or digests are preferred.

The iRibo image includes the compatibility patch and R packages required by
`GenerateTranslatome.R`. Do not substitute the unmodified upstream image.

## Samplesheet creation

```bash
python3 pipelines/translon-analysis/bin/make_samplesheet.py \
  --riboseq-outdir /hps/nobackup/.../riboseq/results \
  --output /hps/nobackup/.../translon/samplesheet.tsv
```

Check the generated sheet before submission. For the full suite it must have
non-empty `transcriptome_bam`, `genome_bam`, `ribo_fastq` and `offsets` columns
for samples that are being sent to the corresponding tools.

## Full-suite launch

Use the `slurm_apptainer` profile when Apptainer is available:

```bash
nextflow run pipelines/translon-analysis \
  -c pipelines/translon-analysis/conf/hpc_apptainer.config \
  -profile slurm_apptainer \
  -resume \
  --samplesheet /hps/nobackup/.../translon/samplesheet.tsv \
  --gtf /hps/nobackup/.../references/genes.gtf \
  --fasta /hps/nobackup/.../references/genome.fa \
  --proteome_fasta /hps/nobackup/.../references/proteome.fa \
  --tools all \
  --ribosomal_fasta /hps/nobackup/.../references/rRNA.fa \
  --adapter_fasta /hps/nobackup/.../references/adapter.fa \
  --orfrater_model /hps/nobackup/.../models/orfrater \
  --outdir /hps/nobackup/.../translon/results
```

The process definitions use immutable public registry image references directly.
Apptainer will cache those images under `NXF_APPTAINER_CACHEDIR`; pre-pulling is
optional but recommended for controlled HPC runs. The custom GHCR images must
exist and be readable from the cluster before launching.

Start with `--tools ribocode,ribotricer,ribotish` on one sample, then add the
legacy and cohort-level callers. PRICE is cohort-level and deliberately uses
more memory; do not use it as the first smoke test.

For the opt-in sharded RiboCode path, add `--partition_count N`
`--partition_mode transcriptome` and, for iRibo genomic sharding, also provide
`--partition_fai /absolute/path/to/genome.fa.fai` with
`--partition_mode genome`. Keep `--partition_padding 0` until equivalence has
been demonstrated on the target reference.

## Acceptance checks

Every tool is accepted only when its Nextflow process exits zero, its native
output is present (or it records a genuine native zero-call result), its
standardised TSV and BED12 exist, and the trace contains the expected image and
resource allocation. Run the adapter tests before submission:

```bash
python3 -m unittest pipelines/translon-analysis/tests/test_standardise_caller.py
git diff --check
```
