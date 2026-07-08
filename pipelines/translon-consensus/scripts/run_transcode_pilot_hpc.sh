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
#   export RUN_CONSENSUS_REPORT=true
#   export RUN_NEXTFLOW=0  # only stage/audit inputs

PILOT_ROOT="${PILOT_ROOT:-/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/full_pilot_results}"
RUN_ROOT="${RUN_ROOT:-/hps/nobackup/flicek/ensembl/genebuild/jackt/riboseq/pilot/translon_consensus_run}"
PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

GENCODE_FASTA="${GENCODE_FASTA:-}"
GENCODE_GTF="${GENCODE_GTF:-}"
GENCODE_FASTA_FAI="${GENCODE_FASTA_FAI:-${GENCODE_FASTA}.fai}"
NF_PROFILE="${NF_PROFILE:-slurm,singularity}"
CONSENSUS_EXCLUDE_TOOLS="${CONSENSUS_EXCLUDE_TOOLS:-RibORF2}"
RUN_CONSENSUS_REPORT="${RUN_CONSENSUS_REPORT:-false}"
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
    local source_feature_class="${6:-unknown}"

    [[ -d "${root}" ]] || return 0
    find "${root}" -maxdepth 1 -type f -name "${pattern}" -print0 |
        sort -z |
        while IFS= read -r -d '' path; do
            name="$(basename "${path}")"
            sample="${name}"
            sample="$(printf '%s' "${sample}" | sed -E "s/${prefix_regex}//; s/${suffix_regex}//")"
            printf '%s\t%s\t%s\t%s\n' "${path}" "${tool}" "${sample}" "${source_feature_class}"
        done
}

{
    printf 'path\ttool\tsample_id\tsource_feature_class\n'

    # Prefer raw annotated/unannotated iRibo outputs so CDS/non-CDS class remains explicit.
    write_rows "iRibo" "${PILOT_ROOT}/iRibo_results" "annotated_orfs_*.bed" '^annotated_orfs_' '\.bed$' "cds"
    write_rows "iRibo" "${PILOT_ROOT}/iRibo_results" "unannotated_orfs_*.bed" '^unannotated_orfs_' '\.bed$' "non_cds"

    # ORFQuant has separate annotated and novel genomic BEDs in bed_files.
    if [[ -d "${PILOT_ROOT}/ORFQuant_results/bed_files" ]]; then
        find "${PILOT_ROOT}/ORFQuant_results/bed_files" -maxdepth 1 -type f \( -name "annotated_orf_exon_genomic_*.bed" -o -name "novel_orf_exon_genomic_*.bed" \) -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                name="$(basename "${path}")"
                class="non_cds"
                [[ "${name}" == annotated_orf_exon_genomic_* ]] && class="cds"
                sample="$(printf '%s' "${name}" |
                    sed -E 's/^(annotated|novel)_orf_exon_genomic_//; s/\.bed$//; s/\.Aligned\.sortedByCoord\.out$//; s/_trimmed\.Aligned\.sortedByCoord\.out$//; s/_S[0-9]+_R1_001_trimmed$//; s/Ribo_pancreas_pooled/Ribo_Pancreas_pooled/')"
                printf '%s\t%s\t%s\t%s\n' "${path}" "ORFQuant" "${sample}" "${class}"
            done
    fi

    # RiboTIE deliverables are split by annotated/novel when present; otherwise fall back to processed ncORFs.
    ribotie_rows=0
    if [[ -d "${PILOT_ROOT}/RiboTIE_results/deliverables" ]]; then
        while IFS= read -r -d '' path; do
            class="unknown"
            case "${path}" in
                */annotated/*) class="cds" ;;
                */novel/*) class="non_cds" ;;
            esac
            name="$(basename "${path}")"
            sample="$(printf '%s' "${name}" |
                sed -E 's/^ncORFs_//; s/^(RiboTIE_Annotations_|RiboTIE-)//; s/\.(bed|bed12|gtf|gff|gff3)$//; s/GENELAB-000/GENELAB/g; s/GENELAB_000/GENELAB/g')"
            printf '%s\t%s\t%s\t%s\n' "${path}" "RiboTIE" "${sample}" "${class}"
            ribotie_rows=$((ribotie_rows + 1))
        done < <(find "${PILOT_ROOT}/RiboTIE_results/deliverables" -type f \( -name "*.bed" -o -name "*.bed12" -o -name "*.gtf" -o -name "*.gff" -o -name "*.gff3" \) -print0 | sort -z)
    fi
    if [[ "${ribotie_rows}" -eq 0 ]]; then
        write_rows "RiboTIE" "${PILOT_ROOT}/Processed outputs/RiboTIE" "*.bed" '^ncORFs_' '\.bed$' "non_cds"
    fi

    # PRICE caller output is split into known and novel files; keep that split for the DB.
    if [[ -d "${PILOT_ROOT}/PRICE_results/price" ]]; then
        write_rows "PRICE" "${PILOT_ROOT}/PRICE_results/price" "*.known.bed" '^' '\.known\.bed$' "cds"
        write_rows "PRICE" "${PILOT_ROOT}/PRICE_results/price" "*.novel.bed" '^' '\.novel\.bed$' "non_cds"
    fi

    # RibORF2 pancreas fastq outputs are split into annotated and novel files.
    if [[ -d "${PILOT_ROOT}/fastq_RibORF2.0_ORFidentification" ]]; then
        find "${PILOT_ROOT}/fastq_RibORF2.0_ORFidentification" \
            -mindepth 2 -maxdepth 2 -type f -name "*.bed" -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                name="$(basename "${path}")"
                class="non_cds"
                [[ "${name}" == *annotatedORFs.bed || "${name}" == *annotatedorfs.bed ]] && class="cds"
                sample="$(printf '%s' "${name}" |
                    sed -E 's/_(annotatedORFs|annotatedorfs|novelsmorf|novelsmorfs)\.bed$//; s/^Pancreas_?([0-9]+)_//; s/^pooled_pancreas$/Ribo_Pancreas_pooled/')"
                printf '%s\t%s\t%s\t%s\n' "${path}" "RibORF2" "${sample}" "${class}"
            done
    fi

    # Include any BED-like RibORF caller files if subdirectories contain them.
    if [[ -d "${PILOT_ROOT}/RibORF_results/RibORF_Output_bamtoORFcalling" ]]; then
        find "${PILOT_ROOT}/RibORF_results/RibORF_Output_bamtoORFcalling" \
            -mindepth 2 -maxdepth 2 -type f \( -name "*.bed" -o -name "*.bed12" \) -print0 |
            sort -z |
            while IFS= read -r -d '' path; do
                sample="$(basename "$(dirname "${path}")" | sed -E 's/_trimmed$//')"
                printf '%s\t%s\t%s\t%s\n' "${path}" "RibORF" "${sample}" "unknown"
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
        source_class = row.get("source_feature_class", "unknown")
        sample_id = normalise_sample(row["sample_id"])
        if source_class in {"cds", "non_cds"}:
            sample_id = f"{sample_id}.{source_class}"
        row["sample_id"] = sample_id
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
import sys
from collections import defaultdict
from pathlib import Path

from translon_db_standardise import normalise_sample

src_root = Path(sys.argv[1])
dest_root = Path(sys.argv[2])
excluded = {part for part in sys.argv[3].split(",") if part}

included = 0
for tool_dir in sorted(src_root.iterdir()):
    if not tool_dir.is_dir() or tool_dir.name in excluded:
        continue
    out_tool = dest_root / tool_dir.name
    out_tool.mkdir(parents=True, exist_ok=True)
    grouped = defaultdict(list)
    for src in sorted(tool_dir.iterdir()):
        sample = normalise_sample(src.stem)
        grouped[sample].append(src)

    for sample, paths in sorted(grouped.items()):
        suffix = paths[0].suffix or ".bed"
        dest = out_tool / f"{sample}{suffix}"
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        with dest.open("wb") as out_handle:
            for src in paths:
                target = src.resolve() if src.is_symlink() else src
                with target.open("rb") as in_handle:
                    out_handle.write(in_handle.read())
                    out_handle.write(b"\n")
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
    --run_consensus_report "${RUN_CONSENSUS_REPORT}" \
    --gencode_fasta "${GENCODE_FASTA}" \
    --gencode_fasta_fai "${GENCODE_FASTA_FAI}" \
    --gencode_gtf "${GENCODE_GTF}" \
    --outdir "${RUN_ROOT}/results" \
    --tracedir "${RUN_ROOT}/results/pipeline_info" \
    -resume \
    2>&1 | tee "${RUN_ROOT}/logs/nextflow.$(date +%Y%m%d_%H%M%S).log"
