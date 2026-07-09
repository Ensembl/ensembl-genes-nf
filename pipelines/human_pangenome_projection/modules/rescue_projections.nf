process RESCUE_PROJECTIONS {
    tag "${meta.id}"
    label 'process_medium'

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
    // Four rescue passes (validation-guided remap, internal realign, splice/codon micro-shift),
    // re-validation, and protein QC. Mutates the mapped gene set.
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hpp rescue-projections \\
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
    touch ${prefix}_state/protein_qc.jsonl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: 0.1.0
    END_VERSIONS
    """
}
