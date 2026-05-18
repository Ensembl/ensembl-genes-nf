# Tool Output Standardisation for the Translon DB

This document records what each supported tool produces natively, what
pre-processing is required before the DB builder can ingest it, and how
`source_feature_class` (CDS vs non-CDS) is assigned.

The canonical entry point for a new annotation set is
`scripts/build_translon_manifest.py`, which produces a three-column TSV
(`path / tool / sample_id`) passed to `translon_db_standardise.py` via
`--manifest` (or `--manifest_tsv` in Nextflow).

---

## How `source_feature_class` is assigned

`translon_db_standardise.py` calls `infer_source_feature_class(path, raw_label)`
on the **original file path** (not a renamed staging copy).  The function
searches the concatenated posix path + stem for keyword patterns:

| Pattern matched | Assigned class |
|---|---|
| `known`, `annotated`, `annotated_orfs`, `cds` | `cds` |
| `novel`, `unannotated`, `ncorfs`, `non_cds` | `non_cds` |
| *(no match)* | `unknown` |

**Critical**: `unannotated` contains the substring `annotated`.  Any code
that classifies by substring must test `unannotated` before `annotated`
(see `_classify()` in `build_translon_manifest.py`).

When the manifest approach is used, original file paths are preserved so
this inference works correctly.  If files are staged/renamed first, the
keywords are lost and everything becomes `unknown`.

---

## Sample name normalisation

All tools go through `normalise()` in `build_translon_manifest.py`, which
applies the following in order:

1. Strip trailing `_1` (e.g. `SRR15513179_1` → `SRR15513179`).
2. Normalise GENELAB accession formatting (`GENELAB-0001026` → `GENELAB1026`).
3. Apply `PRICE_SAMPLE_MAP` — converts PRICE internal names to canonical IDs
   (e.g. `Pancreas_0` → `SRR11005875_to_79`, `Fibo_4` → `Fib_24_45m`).
4. Apply `IRIBO_POOLED_MAP` — expands iRibo abbreviated pooled names
   (e.g. `Ribo_EC_p` → `Ribo_EC_pooled`).
5. Apply `MISC_FIXES` — corrects capitalisation discrepancies across tools
   (e.g. `Ribo_pancreas_pooled` → `Ribo_Pancreas_pooled`).

To add a new annotation set, extend these maps as needed before running
`build_translon_manifest.py`.

---

## Tool-by-tool requirements

### PRICE

| Property | Value |
|---|---|
| Native format | BED12 |
| DB parser | `bed12` |
| CDS/non-CDS split | By filename: `{sample}.known.bed` = CDS, `{sample}.novel.bed` = non-CDS |
| Pre-processing | None — BED12 is read directly |
| Sample name | PRICE internal name in stem (e.g. `Pancreas_0`) → `PRICE_SAMPLE_MAP` |
| Fastq variants | Files with `mymapping` in stem are a FASTQ re-run duplicate — **exclude** |
| Known gaps | `Pancreas_1`, `Endo_1`, `Fibo_1` BED files absent from pilot run (PRICE job failures) |

**Collection pattern:**
```
PRICE_results/price/{sample}.known.bed   → CDS
PRICE_results/price/{sample}.novel.bed   → non-CDS
```
Exclude `*mymapping*.bed`.

---

### RiboTIE

| Property | Value |
|---|---|
| Native format | GTF (CDS features, `ORF_id` attribute) |
| DB parser | `ribotie_gtf` — detected when `source == "RiboTIE"` or `ORF_id` in attributes |
| Coordinate convention | 1-based GFF-style; `converted_intervals` applies `ribotie_gff1_plus_stop` rule (convert to 0-based + extend terminal block by 3 nt to include stop codon) |
| CDS/non-CDS split | By subdirectory and filename: `deliverables/annotated/` = CDS, `deliverables/novel/` = non-CDS |
| Pre-processing | None — GTF read directly |
| Sample name | Strip `db_` prefix, `.annotated.out.gtf` / `.novel.out.gtf` suffix, `.Aligned.toTranscriptome.out` infix, `_Transcriptome` suffix (pooled) |
| Fastq variants | `raw/fastq/db_Pancreas_*.gtf` are from a FASTQ re-run; `deliverables/` Pancreas-named entries (CSV only) are also FASTQ — **exclude both**, use SRR-ID GTF files only |
| Raw CSV files | `RiboTIE_results/raw/db_*.csv` are unprocessed TranslonScorer output — do not use in place of deliverable GTFs |

**Collection pattern:**
```
RiboTIE_results/deliverables/annotated/db_*.annotated.out.gtf   → CDS
RiboTIE_results/deliverables/novel/db_*.novel.out.gtf           → non-CDS
```
Exclude `unfiltered/` subdirectory contents.

---

### ORFQuant

| Property | Value |
|---|---|
| Native format | GFF3 (CDS features, ORF identifier in column 9 — no `=` sign, contains `ENST`) |
| DB parser | `orfquant_gff` — detected when `source == "ORFQuant"` or `ENST` in col 9 without `=` |
| Coordinate convention | `orfquant_bedlike_plus_stop` rule: coordinates already 0-based-like; terminal block extended +3 nt for stop codon |
| CDS/non-CDS split | By filename prefix: `annotated_orf_exon_genomic_` = CDS, `novel_orf_exon_genomic_` = non-CDS |
| Pre-processing | None — GFF3 read directly |
| Sample name | Strip `annotated_orf_exon_genomic_` / `novel_orf_exon_genomic_` prefix, `_S##_R1_001_trimmed` sequencing run suffix, `.Aligned.sortedByCoord.out` suffix |
| Duplicates | Some samples have both `*.bed` and `*.Aligned.sortedByCoord.out.bed` — **prefer the shorter name** (same data, different run-label naming) |
| Pooled name | `Ribo_pancreas_pooled` (lowercase p) → normalised to `Ribo_Pancreas_pooled` via `MISC_FIXES` |

**Collection pattern:**
```
ORFQuant_results/bed_files/annotated_orf_exon_genomic_*.bed   → CDS
ORFQuant_results/bed_files/novel_orf_exon_genomic_*.bed       → non-CDS
```
When two files resolve to the same `(sample_id, class)` key, keep the shorter filename.

---

### iRibo

| Property | Value |
|---|---|
| Native format | GFF3 (CDS features, `ID=candidate_orfNNNNN` attribute) |
| DB parser | `iribo_gff` — detected when `source == "iRibo"` or `ID=` in attributes |
| Coordinate convention | `iribo_bedlike` rule: 1-based GFF coordinates converted to 0-based; **no stop-codon extension** (iRibo does not include the stop) |
| CDS/non-CDS split | By filename prefix: `annotated_orfs_` = CDS, `unannotated_orfs_` = non-CDS |
| Pre-processing | None — GFF3 read directly |
| Sample name | Strip `annotated_orfs_` / `unannotated_orfs_` prefix, `_S##_R1_001` sequencing run suffix |
| Pooled names | Abbreviated: `Ribo_EC_p`, `Ribo_Fib_p`, `Ribo_pancreas_p` → expanded via `IRIBO_POOLED_MAP` |
| Fastq variants | Files ending `_fastq.bed` are from a FASTQ re-run — **exclude** |
| Coverage gap | No individual endothelial samples (`SRR15513179`–`SRR15513182` equivalents absent); only `Ribo_EC_p` (pooled) present |

**Collection pattern:**
```
iRibo_results/annotated_orfs_*.bed     → CDS
iRibo_results/unannotated_orfs_*.bed   → non-CDS
```
Exclude `*_fastq.bed`.

> **Warning**: `unannotated` contains the substring `annotated`.  Classification
> code must test `unannotated` before `annotated` to avoid misclassification.

---

### RibORF2

| Property | Value |
|---|---|
| Native format | GenePred (`.txt`) — **not natively supported by `translon_db_standardise.py`** |
| DB parser | `bed12` (after conversion) |
| CDS/non-CDS split | **None** — RibORF2 outputs a single set of representative validated ORFs (`repre.valid.ORF.genepred.txt`) with no CDS/novel distinction |
| Pre-processing required | Convert GenePred → BED12 with `genePredToBed` (UCSC tools) before building manifest |
| `source_feature_class` | Will be `unknown`; CDS classification can only be done post-hoc via the `reference_cds` JOIN in the DB |
| Sample name | Strip `_S##_R1_001_trimmed` sequencing run suffix from directory name |

**Pre-processing command:**
```bash
RIBORF2_DIR=full_pilot_results/RibORF_results/RibORF_Output_bamtoORFcalling
OUT_DIR=staged_v2/riborf2_converted
mkdir -p $OUT_DIR
for d in $RIBORF2_DIR/*/; do
    sample=$(basename $d | sed 's/_S[0-9]*_R1_001_trimmed//')
    genePredToBed $d/repre.valid.ORF.genepred.txt $OUT_DIR/${sample}.bed12
done
```

**Collection pattern:**
```
riborf2_converted/{sample}.bed12   → class unknown (no keyword in path)
```

---

## Running the full pipeline

### Step 1 — Convert RibORF2 GenePred (once per annotation set)

```bash
bash scripts/convert_riborf2_genepred.sh \
    <riborf2_input_dir> \
    <riborf2_converted_dir>
```

### Step 2 — Build manifest

```bash
python3 scripts/build_translon_manifest.py \
    --results-dir /path/to/full_pilot_results \
    --riborf2-converted /path/to/riborf2_converted \
    --out /path/to/manifest.tsv
```

Review the printed summary.  Expected output for the pilot dataset:

- 25 expected samples (11 Fibroblast, 7 Endothelial, 7 Pancreas including pooled)
- PRICE, RiboTIE, ORFQuant, iRibo: both CDS and non-CDS files per sample
- RibORF2: single file per sample, class `unknown`
- PRICE missing 3 samples (`Pancreas_1`, `Endo_1`, `Fibo_1`) — PRICE job failures

### Step 3 — Run DB build via Nextflow

```bash
nextflow run pipelines/translon-consensus/main.nf \
    -profile slurm,singularity \
    --bed_results_dir   /path/to/staged_v2 \
    --manifest_tsv      /path/to/manifest.tsv \
    --outdir            /path/to/translon_db_output \
    --gencode_fasta     /path/to/GRCh38.fa \
    --gencode_gtf       /path/to/gencode.v45.annotation.gtf.gz \
    --build_translon_db true \
    --run_consensus_report false
```

`bed_results_dir` is still required as a process input path even when
`manifest_tsv` is provided; pass the staging directory (any valid path).

### Step 3 (alternative) — Run DB build standalone

```bash
bsub -M 64000 -q production \
    python3 bin/translon_db_standardise.py \
        --input-root  /path/to/staged_v2 \
        --manifest    /path/to/manifest.tsv \
        --out-dir     /path/to/translon_db_output \
        --fasta       /path/to/GRCh38.fa \
        --gtf         /path/to/gencode.v45.annotation.gtf.gz
```

---

## Adding a new annotation set

1. Check the tool's native output format against the table above.
2. If the format is GenePred or another unsupported type, convert to BED12
   first and name the output files with `known`/`novel`/`annotated`/`unannotated`
   keywords so `source_feature_class` is inferred correctly.
3. Add any new sample name mappings to `PRICE_SAMPLE_MAP`, `IRIBO_POOLED_MAP`,
   or `MISC_FIXES` in `build_translon_manifest.py`.
4. Add a `collect_<tool>()` function following the existing pattern, or extend
   an existing one to handle additional output directories.
5. Add the new expected samples to `EXPECTED_SAMPLES` and re-run
   `build_translon_manifest.py` to validate coverage before the DB build.
