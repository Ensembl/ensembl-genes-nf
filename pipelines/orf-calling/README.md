# ORF Calling Pipeline (Wave 1 + Wave 2)

A Nextflow DSL2 pipeline that runs ORF callers on Ribo-seq preprocessing outputs from `pipelines/riboseq`.

Wave 1 tools (implemented with three-step modules: prep, run, parse):
- RiboCode (transcriptome BAM)
- ribotricer (transcriptome BAM)
- RiboTaper (genome BAM)
- ORFquant (transcriptome BAM)
- Rp-Bp (genome BAM)

Design:
- Each tool module is a single responsibility: prepare inputs, run tool, standardize outputs.
- Python utilities (functional, no classes) in `bin/` for input preparation and output parsing.
- Standardized outputs per tool, per sample: `<sample>.orf_calls.tsv` and `<sample>.orf_calls.bed12`.
- Transcriptome vs genome: callers output in their native space; transcript-space BED12 uses `transcript_id` as chrom (interim), genome-space callers use genomic chrom.

Quick start (stub):

```
nextflow run pipelines/orf-calling/main.nf \
  --manifest pipelines/orf-calling/test_stub/samples.csv \
  --gtf ensembl-genes-nf/test_stub/ann.gtf \
  --fasta ensembl-genes-nf/test_stub/ref.fa \
  --tool ribocode \
  -stub -profile local
```

Artifacts are written under `results/orf_calls/<tool>/`.

Support matrix (Wave 1):
- RiboCode: RUN+PARSE complete; TSV+BED12 emitted; optional length-gating
- ribotricer: RUN+PARSE complete; TSV+BED12 emitted; length-gating via `--read_lengths`
- ORFquant: RUN shim + detection; TSV+BED12 emitted (content depends on tool availability in environment)
- RiboTaper: RUN shim + detection; TSV+BED12 emitted (content depends on tool availability in environment)
- Rp-Bp: RUN shim + detection; TSV+BED12 emitted (content depends on tool availability in environment)

Notes:
- For reproducible runs, prefer `-profile docker` with pinned containers; fall back to conda as needed.

Wave 2 tools (three-step modules with standardized outputs; RUN shims in place):
- iRibo (genome BAM)
- ORF-RATER (transcriptome BAM)
- PRICE (genome BAM)
- RibORF (genome BAM)
- Ribo-TISH (genome BAM)
- RiboTIE (genome BAM)

Usage examples:
- Run all Wave 2 (stub):
  - `nextflow run pipelines/orf-calling/main.nf --manifest pipelines/orf-calling/test_stub/samples.csv --gtf ensembl-genes-nf/test_stub/ann.gtf --fasta ensembl-genes-nf/test_stub/ref.fa --tool all-wave2 -stub -profile local`

Output schema and BED12 expectations are identical to Wave 1.
