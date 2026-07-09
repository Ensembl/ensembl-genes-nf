process FINALISE {
    tag "${meta.id}"
    label 'process_low'

    container "${params.hpp_container}"

    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.{gff3,json,tsv}"

    input:
    tuple val(meta), path(state_in)

    output:
    tuple val(meta), path("${prefix}.mapped.gff3"),       emit: gff
    tuple val(meta), path("${prefix}.stats.json"),        emit: stats
    path "${prefix}.stats.audit.json",   optional: true,  emit: audit
    path "${prefix}.stats.removed.tsv",  optional: true,  emit: removed
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Statistics, synteny analysis, and final GFF3 + stats output generation.
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    hpp finalise \\
        --in-state ${state_in} \\
        --output-gff ${prefix}.mapped.gff3 \\
        --output-stats ${prefix}.stats.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: \$(hpp --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.gff3 ${prefix}.stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: 0.1.0
    END_VERSIONS
    """
}
