process MINIMAP2_ALIGN {
    tag "${meta.id}:alignment"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0'

    input:
    tuple val(meta), path(reads)
    path minimap_index

    output:
    tuple val(meta), path('alignment.sam'), emit: sam
    path 'versions.yml', emit: versions

    script:
    def secondary = task.ext.secondary ?: (params.secondary_mode == 'yes' ? '--secondary=yes' : '--secondary=no')
    def preset = meta.minimap2_preset ?: params.minimap2_preset
    """
    minimap2 -t ${task.cpus} -ax ${preset} ${secondary} "${minimap_index}" "${reads}" > alignment.sam
    test -s alignment.sam || { echo "minimap2 produced no SAM for ${meta.id}" >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    printf '@HD\tVN:1.6\tSO:unsorted\n' > alignment.sam
    printf '"%s":\n    minimap2: 2.28\n' '${task.process}' > versions.yml
    """
}
