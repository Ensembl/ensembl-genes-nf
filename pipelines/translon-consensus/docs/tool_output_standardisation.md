# Tool Output Standardisation for the Translon DB

This document records what each supported tool produces natively, what
pre-processing is required before the DB builder can ingest it, and how the
auditable source manifest is built.

The canonical entry point for a new annotation set is
`scripts/build_translon_manifest.py`, which produces a provenance-rich TSV
passed to `translon_db_standardise.py` via `--manifest` (or `--manifest_tsv`
in Nextflow). Only rows with `ingest_status=matched` are parsed; rows marked
`empty_at_source`, `blocked`, or `excluded` remain in the TSV for audit.

---

## Native class vs central class

Tool-native labels such as PRICE `known/novel`, RiboTIE `annotated/novel`,
iRibo `annotated/unannotated`, and ORFQuant `annotated/novel` are recorded in
the manifest as `native_class`. They are not treated as the central CDS/non-CDS
classification. Manifest rows set `source_feature_class=unknown`, allowing
`translon_db_standardise.py` to assign the central class from the GENCODE
`reference_cds` join.

This matters most for PRICE: GEDI's `d.type.annotated` flag is PRICE's native
annotation call, not the pilot DB's reference-derived classification.

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
| Native class split | By filename: `{sample}.known.bed` = `annotated`, `{sample}.novel.bed` = `novel` |
| Pre-processing | None — BED12 is read directly |
| Sample name | PRICE internal name in stem (e.g. `Pancreas_0`) → `PRICE_SAMPLE_MAP` |
| BAM route | `Pancreas_0..5`, `Pancreas`, `Fibo_*`, and `Endo_*` files |
| FASTQ route | `pancreas_mymapping_0..5` and `pancreas_mymapping` files |
| Indexing | `Pancreas_0` and `pancreas_mymapping_0` map to pilot pancreas replicate 1 (`SRR11005875_to_79`) |

**Collection pattern:**
```
PRICE_results/price/{sample}.known.bed   -> native_class=annotated
PRICE_results/price/{sample}.novel.bed   -> native_class=novel
```

---

### RiboTIE

| Property | Value |
|---|---|
| Native format | CSV deliverables |
| DB parser | `translonscorer_csv` |
| Native class split | By deliverable subdirectory: `annotated/` = `annotated`, `novel/` = `novel` |
| Pre-processing | None — canonical CSV deliverables are read directly |
| Sample name | Strip `db_`, `.annotated.out` / `.novel.out`, `.Aligned.toTranscriptome.out`, and `_Transcriptome` |
| BAM route | SRR/Fib/pooled transcriptome-named deliverables |
| FASTQ route | `db_Pancreas_1..6.*.csv` deliverables |
| Exclusions | GTF sidecars and `unfiltered/` deliverables are retained as excluded rows |

**Collection pattern:**
```
RiboTIE_results/deliverables/annotated/db_*.annotated.out.csv   -> native_class=annotated
RiboTIE_results/deliverables/novel/db_*.novel.out.csv           -> native_class=novel
```

---

### ORFQuant

| Property | Value |
|---|---|
| Native format | GFF3 (CDS features, ORF identifier in column 9 — no `=` sign, contains `ENST`) |
| DB parser | `orfquant_gff` — detected when `source == "ORFQuant"` or `ENST` in col 9 without `=` |
| Coordinate convention | `orfquant_bedlike_plus_stop` rule: coordinates already 0-based-like; terminal block extended +3 nt for stop codon |
| Native class split | By filename prefix: `annotated_orf_exon_genomic_` = `annotated`, `novel_orf_exon_genomic_` = `novel` |
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
| Native class split | By filename prefix: `annotated_orfs_` = `annotated`, `unannotated_orfs_` = `novel` |
| Pre-processing | None — GFF3 read directly |
| Sample name | Strip `annotated_orfs_` / `unannotated_orfs_` prefix, `_S##_R1_001` sequencing run suffix |
| Pooled names | Abbreviated: `Ribo_EC_p`, `Ribo_Fib_p`, `Ribo_pancreas_p` → expanded via `IRIBO_POOLED_MAP` |
| Fastq variants | Files ending `_fastq.bed` are the pancreas FASTQ route and are included with `route=fastq` |

**Collection pattern:**
```
iRibo_results/annotated_orfs_*.bed     -> native_class=annotated
iRibo_results/unannotated_orfs_*.bed   -> native_class=novel
```

> **Warning**: `unannotated` contains the substring `annotated`.  Classification
> code must test `unannotated` before `annotated` to avoid misclassification.

---

### RibORF2

| Property | Value |
|---|---|
| Native format | GenePred (`.txt`) — **not natively supported by `translon_db_standardise.py`** |
| DB parser | `bed12` (after conversion) |
| Native class split | **None** — converted RibORF2 BED12 rows are recorded as `native_class=unknown` |
| Pre-processing required | Convert GenePred → BED12 with `genePredToBed` (UCSC tools) before building manifest |
| `source_feature_class` | `unknown`; central classification is via the `reference_cds` JOIN in the DB |
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
    --out /path/to/translon_manifest.tsv \
    --design-matrix-out /path/to/design_matrix.tsv
```

Review the printed summary and the manifest status columns. The important
status values are:

- `matched`: parsed by `translon_db_standardise.py`.
- `empty_at_source`: expected source route exists conceptually but has no raw
  output, for example RibORF2 FASTQ.
- `blocked`: expected design-matrix cell with no matched source file.
- `excluded`: raw leaf was seen but is outside canonical ingest rules, such as
  RiboTIE unfiltered sidecars or duplicate ORFQuant names.

The manifest builder exits non-zero if the manifest violates the expected
design matrix. It checks that every required tool/sample/route/native_class
cell is represented, no unexpected source cell is missing, matched source paths
exist and have records, matched paths are unique, FASTQ routes are legal, PRICE
zero-indexed pancreas files map to replicate 1, and every known raw leaf is
either matched or excluded with a reason.

For a standalone Stage 1 report, run:

```bash
python3 scripts/validate_translon_manifest.py \
    --manifest /path/to/translon_manifest.tsv \
    --results-dir /path/to/full_pilot_results \
    --riborf2-converted /path/to/riborf2_converted \
    --out-dir /path/to/validation/stage1_manifest
```

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

`translon_db_standardise.py` runs a second datacheck before writing outputs. It
requires every `ingest_status=matched` manifest row to appear once in
`parser_manifest`, every matched parser row to parse with `parser_status=ok`,
sample/route/native_class/parser metadata to match the manifest, and per-file
translon counts in `translons` to agree with `parser_manifest`.

For a standalone Stage 2 report, run:

```bash
python3 scripts/validate_translon_db_against_manifest.py \
    --manifest /path/to/translon_manifest.tsv \
    --db /path/to/translon_db/translons.sqlite \
    --out-dir /path/to/validation/stage2_db_manifest
```

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
   first and preserve any native annotated/novel split as `native_class`.
3. Add any new sample name mappings to `PRICE_SAMPLE_MAP`, `IRIBO_POOLED_MAP`,
   or `MISC_FIXES` in `build_translon_manifest.py`.
4. Add a `collect_<tool>()` function following the existing pattern, or extend
   an existing one to handle additional output directories.
5. Add the new expected samples to `EXPECTED_SAMPLE_METADATA` and re-run
   `build_translon_manifest.py` to validate coverage before the DB build.
