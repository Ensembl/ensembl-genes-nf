#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
PIPELINE="${ROOT}/pipelines/long_read_tama"
EXAMPLE_DIR="${1:-${ROOT}/test/local_example/run}"
PORT="${PORT:-8765}"
mkdir -p "${EXAMPLE_DIR}"

python3 "${PIPELINE}/test/local_example/create_inputs.py" "${EXAMPLE_DIR}" --url "http://127.0.0.1:${PORT}"
python3 -m http.server "${PORT}" --directory "${EXAMPLE_DIR}" >/tmp/long-read-tama-example-http.log 2>&1 &
SERVER_PID=$!
trap 'kill "${SERVER_PID}" 2>/dev/null || true' EXIT

for attempt in {1..20}; do
    curl --fail --silent "http://127.0.0.1:${PORT}/local_ont.fastq.gz" >/dev/null && break
    sleep 0.2
done

cd "${ROOT}"
mkdir -p "${EXAMPLE_DIR}/inventory-only/classification_report"
python3 "${PIPELINE}/bin/normalise_manifest.py" \
    "${EXAMPLE_DIR}/candidate_manifest.csv" \
    "${EXAMPLE_DIR}/inventory-only/normalised_manifest.tsv" \
    "${EXAMPLE_DIR}/inventory-only/manifest_validation.tsv"
python3 "${PIPELINE}/bin/read_input_classification.py" inspect \
    "${EXAMPLE_DIR}/inventory-only/normalised_manifest.tsv" \
    "${EXAMPLE_DIR}/metadata.json" \
    "${EXAMPLE_DIR}/inventory-only/classification_report"

nextflow run pipelines/long_read_tama/main.nf -stub-run -profile stub \
    --approved_manifest "${EXAMPLE_DIR}/approved_manifest.tsv" \
    --fastq_cache_dir "${EXAMPLE_DIR}/fastq-cache" \
    --reference_fasta "${EXAMPLE_DIR}/genome.fa" \
    --shard_mode contig \
    --outdir "${EXAMPLE_DIR}/sharded" \
    -work-dir "${EXAMPLE_DIR}/work-sharded"

nextflow run pipelines/long_read_tama/main.nf -stub-run -profile stub \
    --approved_manifest "${EXAMPLE_DIR}/approved_manifest.tsv" \
    --fastq_cache_dir "${EXAMPLE_DIR}/fastq-cache-none" \
    --reference_fasta "${EXAMPLE_DIR}/genome.fa" \
    --shard_mode none \
    --outdir "${EXAMPLE_DIR}/unsharded" \
    -work-dir "${EXAMPLE_DIR}/work-unsharded"

if command -v minimap2 >/dev/null && command -v samtools >/dev/null; then
    mkdir -p "${EXAMPLE_DIR}/real-sharding/shards"
    minimap2 -a -x splice "${EXAMPLE_DIR}/genome.fa" "${EXAMPLE_DIR}/local_ont.fastq.gz" \
        | samtools sort -o "${EXAMPLE_DIR}/real-sharding/aligned.sorted.bam"
    samtools index "${EXAMPLE_DIR}/real-sharding/aligned.sorted.bam"
    python3 "${PIPELINE}/bin/inspect_bam_workload.py" \
        "${EXAMPLE_DIR}/real-sharding/aligned.sorted.bam" \
        "${EXAMPLE_DIR}/real-sharding/aligned.sorted.bam.bai" \
        "${EXAMPLE_DIR}/real-sharding/contig_workload.tsv" 1
    python3 "${PIPELINE}/bin/split_bam_by_contig.py" \
        "${EXAMPLE_DIR}/real-sharding/aligned.sorted.bam" \
        "${EXAMPLE_DIR}/real-sharding/contig_workload.tsv" \
        SRR900001 "${EXAMPLE_DIR}/real-sharding/shards" \
        "${EXAMPLE_DIR}/real-sharding/contig_manifest.tsv"
    echo "Real minimap2/samtools sharding smoke test completed"
else
    echo "Skipping real sharding smoke test (minimap2 and samtools are required)"
fi

echo "Local example completed. Inputs and reports are in ${EXAMPLE_DIR}"
