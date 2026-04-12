#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")"/../../.. && pwd)"
PIPE="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/main.nf"
CFG="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/nextflow.config"
MAN="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/test_stub/samples.csv"
GTF="${ROOT}/ensembl-genes-nf/test_stub/ann.gtf"
FA="${ROOT}/ensembl-genes-nf/test_stub/ref.fa"
OUT="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/test_stub/results_ci"

rm -rf "$OUT"

echo "[NF] Wave 1 (stub)"
nextflow run "$PIPE" -c "$CFG" -profile local -stub \
  --manifest "$MAN" --gtf "$GTF" --fasta "$FA" --tool all-wave1 --outdir "$OUT/wave1"

echo "[CHK] Wave 1 artifacts"
grep -q "\torf_calls\.tsv$" <(find "$OUT/wave1" -name '*.orf_calls.tsv')
grep -q "\torf_calls\.bed12$" <(find "$OUT/wave1" -name '*.orf_calls.bed12')

echo "[NF] Wave 2 (stub)"
nextflow run "$PIPE" -c "$CFG" -profile local -stub \
  --manifest "$MAN" --gtf "$GTF" --fasta "$FA" --tool all-wave2 --outdir "$OUT/wave2"

echo "[CHK] Wave 2 artifacts"
grep -q "\torf_calls\.tsv$" <(find "$OUT/wave2" -name '*.orf_calls.tsv')
grep -q "\torf_calls\.bed12$" <(find "$OUT/wave2" -name '*.orf_calls.bed12')

echo "OK: stub tests passed"

