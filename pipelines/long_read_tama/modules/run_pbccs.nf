process RUN_PBCCS {
    tag "${meta.id}:ccs"
    label 'process_high_memory'
    container { params.ccs_container ?: 'https://depot.galaxyproject.org/singularity/pbccs:6.4.0--h9ee0642_0' }

    input:
    tuple val(meta), path(input_artifacts)

    output:
    tuple val(meta), path('consensus.bam'), emit: bam
    path 'versions.yml', emit: versions

    script:
    def ccs_args = task.ext.ccs_args ?: params.ccs_args ?: ''
    """
    input_bam=\$(find . -maxdepth 1 -type f -name '*.bam' -print -quit)
    test -n "\${input_bam}" || { echo "Expected raw-subread BAM for ${meta.id}" >&2; exit 1; }
    ${params.ccs_command ?: 'ccs'} ${ccs_args} "\${input_bam}" consensus.bam
    test -s consensus.bam || { echo "CCS produced no BAM for ${meta.id}" >&2; exit 1; }
    ccs_version=\$(${params.ccs_command ?: 'ccs'} --version 2>&1 | head -n1)
    if [ -n "${params.ccs_expected_version ?: ''}" ]; then
        echo "\${ccs_version}" | grep -F "${params.ccs_expected_version}" >/dev/null || { echo "Unexpected CCS version for ${meta.id}: \${ccs_version}" >&2; exit 1; }
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ccs: \${ccs_version}
    END_VERSIONS
    """

    stub:
    """
    touch consensus.bam
    printf '"%s":\n    ccs: stub\n' '${task.process}' > versions.yml
    """
}
