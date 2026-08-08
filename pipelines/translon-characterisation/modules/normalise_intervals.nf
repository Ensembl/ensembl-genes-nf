process NORMALISE_INTERVALS {
    label 'process_low'
    container params.ucsc_container
    tag "${meta.id}"
    publishDir "${params.outdir}/00_inputs", mode: 'copy', pattern: '*.bed'

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path('*.bed'), emit: intervals
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def lower = intervals.name.toLowerCase()
    if (lower.endsWith('.bb') || lower.endsWith('.bigbed')) {
        """
        bigBedToBed ${intervals} ${prefix}.bed
        test -s ${prefix}.bed
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ucsc_bigBedToBed: \$(bigBedToBed 2>&1 | head -1 | sed 's/^bigBedToBed v//')
        END_VERSIONS
        """
    } else {
        """
        cp ${intervals} ${prefix}.bed
        test -s ${prefix}.bed
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ucsc_bigBedToBed: not-required
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc_bigBedToBed: 482
    END_VERSIONS
    """
}
