process MINIMAP2_ALIGN {
    tag "${meta.id}:alignment"
    label 'process_high'
    // Keep minimap2 and samtools in one image so SAM streams directly into
    // samtools sort instead of being materialized as a large work-dir file.
    container 'docker://niemasd/minimap2_samtools:2.28_1.20'

    input:
    tuple val(meta), path(reads)
    path minimap_index

    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: bam
    path 'versions.yml', emit: versions

    script:
    def secondary = task.ext.secondary ?: (params.secondary_mode == 'yes' ? '--secondary=yes' : '--secondary=no')
    def preset = meta.minimap2_preset ?: params.minimap2_preset
    """
    minimap2 -t ${task.cpus} -ax ${preset} ${secondary} "${minimap_index}" "${reads}" \\
        | samtools sort -@ ${task.cpus} -m ${params.samtools_sort_memory} -O BAM -o "${meta.id}.sorted.bam" -
    test -s "${meta.id}.sorted.bam" || { echo "alignment produced no BAM for ${meta.id}" >&2; exit 1; }
    samtools index "${meta.id}.sorted.bam"
    test -s "${meta.id}.sorted.bam.bai" || { echo "alignment produced no BAI for ${meta.id}" >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    printf '"%s":\n    minimap2: 2.28\n    samtools: 1.20\n' '${task.process}' > versions.yml
    """
}
