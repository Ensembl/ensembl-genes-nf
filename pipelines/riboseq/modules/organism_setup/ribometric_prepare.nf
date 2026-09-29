/*
 * RIBOMETRIC_PREPARE
 * Prepare annotation file for RiboMetric
 */

process RIBOMETRIC_PREPARE {
    tag "${organism}_${version}"
    label 'process_medium'

    container "ghcr.io/lapti-ucc/riboseqorg-nf-ribometric:latest"

    input:
    path(gtf)
    path(fasta)
    val(organism)
    val(version)

    output:
    path "*.tsv",          emit: ribometric_tsv
    path "versions.yml",   emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    RiboMetric prepare -p ${task.cpus} -g ${gtf} ${args}

    # Rename output to include organism/version
    for f in *.tsv; do
        if [ "\$f" != "${prefix}_ribometric.tsv" ]; then
            mv "\$f" "${prefix}_ribometric.tsv" 2>/dev/null || true
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: \$(RiboMetric --version 2>&1 | sed 's/RiboMetric //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${organism}_${version}"
    """
    touch ${prefix}_ribometric.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribometric: 0.1.0
    END_VERSIONS
    """
}
