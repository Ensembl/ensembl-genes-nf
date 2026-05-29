#!/usr/bin/env bash
# Rebuild the translon DB directly from raw pilot results (no Nextflow staging).
#
# Usage:
#   export GENCODE_FASTA=/path/to/GRCh38.fa
#   export GENCODE_GTF=/path/to/gencode.v44.annotation.gtf.gz
#   bash rebuild_translon_db.sh
#
# Optional overrides:
#   export PILOT_ROOT=...   # default: /hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/full_pilot_results
#   export RIBORF2_CONVERTED=...  # dir of RibORF2 .bed12 files (from genePredToBed)
#   export OUT_DIR=...      # default: ./translon_db_rebuild

set -euo pipefail

PILOT_ROOT="${PILOT_ROOT:-/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/full_pilot_results}"
RIBORF2_CONVERTED="${RIBORF2_CONVERTED:-${PILOT_ROOT}/riborf2_converted_bed12}"
OUT_DIR="${OUT_DIR:-$(pwd)/translon_db_rebuild}"
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${PIPELINE_DIR}/bin"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

if [[ -z "${GENCODE_FASTA:-}" ]]; then
    echo "ERROR: Set GENCODE_FASTA to the genome FASTA." >&2
    exit 2
fi
if [[ -z "${GENCODE_GTF:-}" ]]; then
    echo "ERROR: Set GENCODE_GTF to the annotation GTF." >&2
    exit 2
fi

echo "=== Translon DB rebuild ==="
echo "PILOT_ROOT   : ${PILOT_ROOT}"
echo "OUT_DIR      : ${OUT_DIR}"
echo "GENCODE_FASTA: ${GENCODE_FASTA}"
echo "GENCODE_GTF  : ${GENCODE_GTF}"
echo

mkdir -p "${OUT_DIR}"
MANIFEST="${OUT_DIR}/translon_manifest.tsv"
DESIGN_MATRIX="${OUT_DIR}/design_matrix.tsv"
VALIDATION_DIR="${OUT_DIR}/validation"

# Step 1: build manifest from original paths
echo "--- Step 1: generating manifest ---"
python3 "${SCRIPTS_DIR}/build_translon_manifest.py" \
    --results-dir "${PILOT_ROOT}" \
    --riborf2-converted "${RIBORF2_CONVERTED}" \
    --out "${MANIFEST}" \
    --design-matrix-out "${DESIGN_MATRIX}"
echo "Manifest datachecks passed."

python3 "${SCRIPTS_DIR}/validate_translon_manifest.py" \
    --manifest "${MANIFEST}" \
    --results-dir "${PILOT_ROOT}" \
    --riborf2-converted "${RIBORF2_CONVERTED}" \
    --out-dir "${VALIDATION_DIR}/stage1_manifest" \
    --ignore-prefix "${OUT_DIR}"
echo "Standalone Stage 1 manifest validation passed."

echo
echo "--- Step 2: building translon DB and checking DB against manifest ---"
PYTHONPATH="${BIN_DIR}${PYTHONPATH:+:${PYTHONPATH}}" python3 "${BIN_DIR}/translon_db_standardise.py" \
    --input-root "${PILOT_ROOT}" \
    --manifest "${MANIFEST}" \
    --out-dir "${OUT_DIR}/translon_db" \
    --fasta "${GENCODE_FASTA}" \
    --gtf "${GENCODE_GTF}"
echo "DB/manifest datachecks passed."

python3 "${SCRIPTS_DIR}/validate_translon_db_against_manifest.py" \
    --manifest "${MANIFEST}" \
    --db "${OUT_DIR}/translon_db/translons.sqlite" \
    --out-dir "${VALIDATION_DIR}/stage2_db_manifest"
echo "Standalone Stage 2 DB/manifest validation passed."

echo
echo "=== Done ==="
echo "DB files: ${OUT_DIR}/translon_db/"
echo "  ${MANIFEST}"
echo "  ${DESIGN_MATRIX}"
echo "  ${VALIDATION_DIR}/"
echo "  translons.sqlite"
echo "  translons.tsv.gz, translon_blocks.tsv.gz"
echo "  parser_manifest.tsv.gz (check parser_status for any unsupported files)"
echo "  cds_recall_by_tool.tsv.gz"
echo
# Quick sanity check
python3 - "${OUT_DIR}/translon_db/translons.sqlite" <<'PY'
import sqlite3, sys
con = sqlite3.connect(sys.argv[1])
rows = con.execute(
    "SELECT source_tool, source_feature_class, COUNT(*) FROM translons "
    "GROUP BY source_tool, source_feature_class ORDER BY source_tool, source_feature_class"
).fetchall()
print(f"{'Tool':<12} {'Class':<10} {'Translons':>10}")
print("-" * 35)
for tool, cls, n in rows:
    print(f"{tool:<12} {cls:<10} {n:>10,}")
print()
bad = con.execute(
    "SELECT detected_tool, COUNT(*) FROM parser_manifest WHERE parser_status!='ok' GROUP BY detected_tool"
).fetchall()
if bad:
    print("WARNING: unsupported/error files:")
    for tool, n in bad:
        print(f"  {tool}: {n} files")
PY
