# European hedgehog Kinnex minimap2/TAMA comparison

This run compares the two Kinnex full-length transcript datasets associated
with specimen `SAN00002332`:

| run | tissue | ENA accession | BioSample | Ensembl assembly |
| --- | --- | --- | --- | --- |
| liver | liver | `ERR13757840` | `SAMEA13207749` | `HEDGEHOG` |
| brain | brain | `ERR13757843` | `SAMEA13207752` | `HEDGEHOG` |

Both runs belong to study `PRJEB60668`, have PacBio Revio metadata, and use
the same Ensembl assembly (`GCA_950295315.1`, Ensembl assembly name
`HEDGEHOG`). The FASTQ files are the ENA-normalised processed products and
their headers are archive-renamed, so the approved manifest uses the
`PACBIO_PROCESSED_FASTQ` route with `splice:hq`.

The reference manifest records the Ensembl release-116 genome and annotation
URLs. Use the unmasked top-level FASTA for minimap2/TAMA; retain the GTF for
downstream reference-guided or annotation-QC work. Do not substitute the
NCBI FASTA if the comparison is intended to remain an Ensembl-reference run.

## Prepare the reference

Download the files listed in `reference_manifest.tsv` into a persistent
directory, decompress the genome FASTA, and build its `.fai` and minimap2
index with the same minimap2/samtools versions used by the pipeline. For
example:

```bash
REF_DIR=/shared/references/erinaceus_europaeus/HEDGEHOG/ensembl-116
mkdir -p "$REF_DIR"
awk -F '\t' 'NR > 1 && $3 == "download" {print $4 "\t" $5}' reference_manifest.tsv |
while IFS=$'\t' read -r url name; do
  curl -fL --retry 3 -o "$REF_DIR/$name" "$url"
done
gunzip -k "$REF_DIR/Erinaceus_europaeus.HEDGEHOG.dna.toplevel.fa.gz"
samtools faidx "$REF_DIR/Erinaceus_europaeus.HEDGEHOG.dna.toplevel.fa"
minimap2 -d "$REF_DIR/Erinaceus_europaeus.HEDGEHOG.dna.toplevel.mmi" \
  "$REF_DIR/Erinaceus_europaeus.HEDGEHOG.dna.toplevel.fa"
```

Check the downloaded files against the release `CHECKSUMS` files before
starting production. The FASTA and GTF URLs are intentionally pinned to
release 116 rather than `current_*` aliases.

## Run

Set `FASTQ_CACHE` to a shared persistent cache. The first run downloads two
FASTQs (approximately the ENA sizes for these accessions) and then reuses
checksum-validated copies on restart.

```bash
FASTQ_CACHE=/shared/caches/ena-fastq
REF_DIR=/shared/references/erinaceus_europaeus/HEDGEHOG/ensembl-116

nextflow run pipelines/long_read_tama/main.nf \
  -profile slurm \
  --approved_manifest pipelines/long_read_tama/runs/erinaceus_europaeus_kinnex/approved_run_manifest.tsv \
  --fastq_cache_dir "$FASTQ_CACHE" \
  --reference_fasta "$REF_DIR/Erinaceus_europaeus.HEDGEHOG.dna.toplevel.fa" \
  --cohort_id erinaceus_europaeus_kinnex \
  --model_backend tama \
  --backend_merge_mode native \
  --minimap2_preset splice:hq \
  --secondary_mode no \
  --shard_mode contig \
  --outdir results/erinaceus_europaeus_kinnex
```

This is a TAMA-only comparison run over a shared minimap2 alignment policy.
The two rows are related tissues from one specimen, not independent specimen
replicates; interpret tissue differences accordingly.

