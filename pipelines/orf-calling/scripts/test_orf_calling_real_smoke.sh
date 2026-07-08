#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   MANIFEST=/abs/path/samples_real.csv GTF=/abs/path/anno.gtf FASTA=/abs/path/ref.fa \
#     bash ensembl-genes-nf/pipelines/orf-calling/scripts/test_orf_calling_real_smoke.sh

ROOT="$(cd "$(dirname "$0")"/../../.. && pwd)"
PIPE="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/main.nf"
CFG="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/nextflow.config"
OUT="${ROOT}/ensembl-genes-nf/pipelines/orf-calling/test_real"

: "${MANIFEST:?Set MANIFEST to a CSV with sample_id,bam,bam_type}"
: "${GTF:?Set GTF to a GTF file}"
: "${FASTA:?Set FASTA to a FASTA file}"

mkdir -p "$OUT"

for TOOL in ribocode ribotricer orfquant ribotaper rpbp; do
  echo "[NF] Real smoke: $TOOL"
  nextflow run "$PIPE" -c "$CFG" -profile conda,local \
    --manifest "$MANIFEST" --gtf "$GTF" --fasta "$FASTA" --tool "$TOOL" --outdir "$OUT/$TOOL"
  echo "[CHK] $TOOL artifacts"
  find "$OUT/$TOOL" -name '*.orf_calls.tsv' -size +0c | grep -q .
  find "$OUT/$TOOL" -name '*.orf_calls.bed12' -size +0c | grep -q .
done

echo "OK: real smoke tests completed"

