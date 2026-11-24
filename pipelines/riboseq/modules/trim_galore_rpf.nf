/*
 * TRIM_GALORE_RPF - Traditional single-end Ribo-Seq adapter trimming.
 *
 * The adapter file is paired with each sample.  Trim Galore 2.x accepts
 * repeated --adapter arguments, so all records in a supplied FASTA are
 * preserved.
 */

process TRIM_GALORE_RPF {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::trim-galore=2.6.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f661248fb4305d2f4481faed98af6fff222a2ee73b23fe8518244d82aaf151e6/data' :
        'community.wave.seqera.io/library/trim-galore:2.6.0--da99b872b7cff766' }"

    input:
    tuple val(meta), path(reads), path(adapter_fasta)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: trimmed_fastq
    tuple val(meta), path("*_trim_galore.log"), emit: log
    tuple val(meta), path("*_trimming_report.txt"), emit: report_text, optional: true
    tuple val(meta), path("*_trimming_report.json"), emit: report_json, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_len = params.min_read_length ?: 20
    def max_len = params.max_read_length ?: 0
    def max_len_arg = max_len > 0 ? "--max_length ${max_len}" : ''
    def input_name = reads.name
    def quality = params.traditional_trim_quality ?: 0

    """
    set -euo pipefail

    adapter_args=''
    if [ "${adapter_fasta.name}" != "NO_ADAPTER_FILE" ]; then
        adapter_args=\$(awk '
            /^>/ {
                if (seq != "") print "--adapter " seq
                seq=""
                next
            }
            {
                gsub(/[[:space:]]/, "")
                seq=seq \$0
            }
            END {
                if (seq != "") print "--adapter " seq
            }
        ' "${adapter_fasta}" | paste -sd' ' -)
        if [ -z "\$adapter_args" ]; then
            echo "Adapter FASTA ${adapter_fasta} contains no sequence." >&2
            exit 2
        fi
    elif [ -n "${params.fastp_adapter_sequence ?: ''}" ]; then
        adapter_args="--adapter ${params.fastp_adapter_sequence}"
    fi

    trim_galore \\
        --length ${min_len} \\
        --quality ${quality} \\
        ${max_len_arg} \\
        \$adapter_args \\
        --gzip \\
        --output_dir . \\
        ${args} \\
        ${reads} > ${prefix}_trim_galore.log 2>&1

    trimmed_output=\$(find . -maxdepth 1 -type f \( -name '*_trimmed.fq.gz' -o -name '*_trimmed.fastq.gz' \) -print -quit)
    if [ -z "\$trimmed_output" ]; then
        echo "Trim Galore did not produce a trimmed FASTQ for ${input_name}." >&2
        exit 3
    fi
    mv "\$trimmed_output" ${prefix}_trimmed.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trim_galore: \$(trim_galore --version | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trimmed.fastq.gz
    touch ${prefix}_trim_galore.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trim_galore: 2.6.0
    END_VERSIONS
    """
}
