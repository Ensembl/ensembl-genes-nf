#!/usr/bin/env bash
set -euo pipefail

# Stage the messy TransCODE pilot caller outputs and run translon-consensus.
#
# Required on HPC:
#   export GENCODE_FASTA=/path/to/genome.fa
#   export GENCODE_GTF=/path/to/annotation.gtf.gz
# Optional:
#   export PILOT_ROOT=/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/full_pilot_results
#   export RUN_ROOT=/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/translon_consensus_run
#   export PIPELINE_DIR=/path/to/ensembl-genes-nf/pipelines/translon-consensus
#   export NF_PROFILE=slurm,singularity
#   export CONSENSUS_EXCLUDE_TOOLS=RibORF2
#   export RUN_NEXTFLOW=0  # only stage/audit inputs

PILOT_ROOT="${PILOT_ROOT:-/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/full_pilot_results}"
RUN_ROOT="${RUN_ROOT:-/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/translon_consensus_run}"
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

GENCODE_FASTA="${GENCODE_FASTA:-}"
GENCODE_GTF="${GENCODE_GTF:-}"
GENCODE_FASTA_FAI="${GENCODE_FASTA_FAI:-${GENCODE_FASTA}.fai}"
NF_PROFILE="${NF_PROFILE:-slurm,singularity}"
CONSENSUS_EXCLUDE_TOOLS="${CONSENSUS_EXCLUDE_TOOLS:-RibORF2}"
RUN_NEXTFLOW="${RUN_NEXTFLOW:-1}"

if [[ ! -d "${PILOT_ROOT}" ]]; then
    echo "PILOT_ROOT does not exist: ${PILOT_ROOT}" >&2
    exit 2
fi

if [[ ! -d "${PIPELINE_DIR}" ]]; then
    echo "PIPELINE_DIR does not exist: ${PIPELINE_DIR}" >&2
    exit 2
fi

if [[ -z "${GENCODE_FASTA}" ]]; then
    echo "Set GENCODE_FASTA to the genome FASTA used for the pilot coordinates." >&2
    exit 2
fi

if [[ -z "${GENCODE_GTF}" ]]; then
    echo "Set GENCODE_GTF to the GENCODE/Ensembl GTF used for annotation." >&2
    exit 2
fi

if [[ ! -f "${GENCODE_FASTA}" ]]; then
    echo "GENCODE_FASTA does not exist: ${GENCODE_FASTA}" >&2
    exit 2
fi

if [[ ! -f "${GENCODE_FASTA_FAI}" ]]; then
    echo "GENCODE_FASTA_FAI does not exist: ${GENCODE_FASTA_FAI}" >&2
    echo "Create it with: samtools faidx ${GENCODE_FASTA}" >&2
    exit 2
fi

if [[ ! -f "${GENCODE_GTF}" ]]; then
    echo "GENCODE_GTF does not exist: ${GENCODE_GTF}" >&2
    exit 2
fi

mkdir -p "${RUN_ROOT}"/{manifests,logs,results,work}

FILE_LIST_RAW="${RUN_ROOT}/manifests/transcode_pilot_orf_inputs.raw.tsv"
FILE_LIST="${RUN_ROOT}/manifests/transcode_pilot_orf_inputs.tsv"
STAGE_DIR="${RUN_ROOT}/staged"
CONSENSUS_INPUT_DIR="${RUN_ROOT}/staged/consensus_orf_inputs"
MERGED_INPUTS="${RUN_ROOT}/manifests/merged_inputs"
rm -rf "${MERGED_INPUTS}"
mkdir -p "${MERGED_INPUTS}"

write_rows() {
    local tool="$1"
    local root="$2"
    local pattern="$3"
    local prefix_regex="$4"
    local suffix_regex="$5"

    [[ -d "${root}" ]] || return 0
    find "${root}" -maxdepth 1 -type f -name "${pattern}" -print0 |
        sort -z |
        while IFS= read -r -d '' path; do
            name="$(basename "${path}")"
            sample="${name}"
            sample="$(printf '%s' "${sample}" | sed -E "s/${prefix_regex}//; s/${suffix_regex}//")"
            printf '%s\t%s\t%s\n' "${path}" "${tool}" "${sample}"
        done
}

{
    printf 'path\ttool\tsample_id\n'

    # Curated processed BED outputs already produced for the pilot.
    write_rows "iRibo" "${PILOT_ROOT}/Processed outputs/iRibo" "*.bed" '^ncORFs_' '\.bed$'
    write_rows "ORFQuant" "${PILOT_ROOT}/Processed outputs/ORFQuant" "*.bed" '^ncORFs_' '\.bed$'
    write_rows "RiboTIE" "${PILOT_ROOT}/Processed outputs/RiboTIE" "*.bed" '^ncORFs_' '\.bed$'

    # PRICE processed directory is empty in the current pilot tree, so merge caller
    # known+novel split files per sample before staging.
    PRICE_MERGED="${MERGED_INPUTS}/PRICE"
    mkdir -p "${PRICE_MERGED}"
    if [[ -d "${PILOT_ROOT}/PRICE_results/price" ]]; then
        find "${PILOT_ROOT}/PRICE_results/price" -maxdepth 1 -type f \( -name "*.known.bed" -o -name "*.novel.bed" \) -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                name="$(basename "${path}")"
                sample="$(printf '%s' "${name}" | sed -E 's/\.(known|novel)\.bed$//')"
                cat "${path}" >> "${PRICE_MERGED}/${sample}.bed"
            done
        write_rows "PRICE" "${PRICE_MERGED}" "*.bed" '^' '\.bed$'
    fi

    # RibORF2 pancreas fastq outputs are useful but live outside Processed outputs.
    RIBORF2_MERGED="${MERGED_INPUTS}/RibORF2"
    mkdir -p "${RIBORF2_MERGED}"
    if [[ -d "${PILOT_ROOT}/fastq_RibORF2.0_ORFidentification" ]]; then
        find "${PILOT_ROOT}/fastq_RibORF2.0_ORFidentification" \
            -mindepth 2 -maxdepth 2 -type f -name "*.bed" -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                name="$(basename "${path}")"
                sample="$(printf '%s' "${name}" |
                    sed -E 's/_(annotatedORFs|annotatedorfs|novelsmorf|novelsmorfs)\.bed$//; s/^Pancreas_?([0-9]+)_//; s/^pooled_pancreas$/Ribo_Pancreas_pooled/')"
                cat "${path}" >> "${RIBORF2_MERGED}/${sample}.bed"
            done
        write_rows "RibORF2" "${RIBORF2_MERGED}" "*.bed" '^' '\.bed$'
    fi

    # Include any BED-like RibORF caller files if subdirectories contain them.
    if [[ -d "${PILOT_ROOT}/RibORF_results/RibORF_Output_bamtoORFcalling" ]]; then
        find "${PILOT_ROOT}/RibORF_results/RibORF_Output_bamtoORFcalling" \
            -mindepth 2 -maxdepth 2 -type f \( -name "*.bed" -o -name "*.bed12" \) -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                sample="$(basename "$(dirname "${path}")" | sed -E 's/_trimmed$//')"
                printf '%s\t%s\t%s\n' "${path}" "RibORF" "${sample}"
            done
    fi
} > "${FILE_LIST_RAW}"

PYTHONPATH="${PIPELINE_DIR}/bin${PYTHONPATH:+:${PYTHONPATH}}" python3 - "${FILE_LIST_RAW}" "${FILE_LIST}" <<'PY'
import csv
import sys

from translon_db_standardise import normalise_sample

src, dest = sys.argv[1:3]
with open(src, newline="") as in_handle, open(dest, "w", newline="") as out_handle:
    reader = csv.DictReader(in_handle, delimiter="\t")
    writer = csv.DictWriter(out_handle, delimiter="\t", fieldnames=reader.fieldnames, lineterminator="\n")
    writer.writeheader()
    for row in reader:
        row["sample_id"] = normalise_sample(row["sample_id"])
        writer.writerow(row)
PY

echo "Wrote file list: ${FILE_LIST}"
echo "Raw file list before sample normalisation: ${FILE_LIST_RAW}"
echo "Input rows by tool:"
awk -F '\t' 'NR > 1 { count[$2]++ } END { for (tool in count) print tool, count[tool] }' "${FILE_LIST}" | sort

rm -rf "${STAGE_DIR}"
python3 "${PIPELINE_DIR}/bin/prepare_orf_result_tree.py" \
    --file-list "${FILE_LIST}" \
    --out-dir "${STAGE_DIR}"

echo "Staging audit:"
awk -F '\t' '
    NR == 1 { next }
    { sub(/\r$/, "", $8); by_status[$8]++; by_tool[$3 "\t" $8]++ }
    END {
        print "stage_status counts:";
        for (status in by_status) print "  " status, by_status[status];
        print "tool/status counts:";
        for (key in by_tool) print "  " key, by_tool[key];
    }
' "${STAGE_DIR}/clean_orf_inputs_manifest.tsv"

echo "Staged tree: ${STAGE_DIR}/clean_orf_inputs"
echo "Detailed manifest: ${STAGE_DIR}/clean_orf_inputs_manifest.tsv"

rm -rf "${CONSENSUS_INPUT_DIR}"
mkdir -p "${CONSENSUS_INPUT_DIR}"
PYTHONPATH="${PIPELINE_DIR}/bin${PYTHONPATH:+:${PYTHONPATH}}" python3 - \
    "${STAGE_DIR}/clean_orf_inputs" \
    "${CONSENSUS_INPUT_DIR}" \
    "${CONSENSUS_EXCLUDE_TOOLS}" <<'PY'
import os
import shutil
import sys
from pathlib import Path

src_root = Path(sys.argv[1])
dest_root = Path(sys.argv[2])
excluded = {part for part in sys.argv[3].split(",") if part}

included = 0
for tool_dir in sorted(src_root.iterdir()):
    if not tool_dir.is_dir() or tool_dir.name in excluded:
        continue
    out_tool = dest_root / tool_dir.name
    out_tool.mkdir(parents=True, exist_ok=True)
    for src in sorted(tool_dir.iterdir()):
        dest = out_tool / src.name
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        target = src.resolve() if src.is_symlink() else src
        os.symlink(target, dest)
        included += 1

print(f"Wrote consensus input tree: {dest_root}")
print(f"Excluded consensus tools: {','.join(sorted(excluded)) or '<none>'}")
print(f"Consensus input files: {included}")
PY

if [[ "${RUN_NEXTFLOW}" == "0" ]]; then
    echo "RUN_NEXTFLOW=0, stopping after staging."
    exit 0
fi

cd "${PIPELINE_DIR}"

nextflow run main.nf \
    -profile "${NF_PROFILE}" \
    -work-dir "${RUN_ROOT}/work" \
    --bed_results_dir "${STAGE_DIR}/clean_orf_inputs" \
    --consensus_bed_results_dir "${CONSENSUS_INPUT_DIR}" \
    --gencode_fasta "${GENCODE_FASTA}" \
    --gencode_fasta_fai "${GENCODE_FASTA_FAI}" \
    --gencode_gtf "${GENCODE_GTF}" \
    --outdir "${RUN_ROOT}/results" \
    --tracedir "${RUN_ROOT}/results/pipeline_info" \
    -resume \
    2>&1 | tee "${RUN_ROOT}/logs/nextflow.$(date +%Y%m%d_%H%M%S).log"
