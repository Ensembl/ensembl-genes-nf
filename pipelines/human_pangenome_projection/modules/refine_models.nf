process REFINE_MODELS {
    tag "${meta.id}"
    label 'process_low'

    container "${params.hpp_container}"

    input:
    tuple val(meta), path(state_in), path(target_fasta)
    path reference

    output:
    tuple val(meta), path("${prefix}_state"), emit: state
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Conflict resolution / refinement (CNV selection, paralog conflicts) and per-gene audit traces.
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hpp refine-models \\
        --in-state ${state_in} \\
        --out-state ${prefix}_state \\
        --ref-fasta ${reference} \\
        --target-fasta ${target_fasta} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: \$(hpp --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_state
    touch ${prefix}_state/refinement_report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: 0.1.0
    END_VERSIONS
    """
}
