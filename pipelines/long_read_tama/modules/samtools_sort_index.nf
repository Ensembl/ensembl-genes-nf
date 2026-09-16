process SAMTOOLS_SORT_INDEX {
    tag "${meta.id}:sort-index"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: bam
    path 'versions.yml', emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -m ${params.samtools_sort_memory} -O BAM -o "${meta.id}.sorted.bam" "${sam}"
    test -s "${meta.id}.sorted.bam" || { echo "samtools sort produced no BAM for ${meta.id}" >&2; exit 1; }
    samtools index "${meta.id}.sorted.bam"
    test -s "${meta.id}.sorted.bam.bai" || { echo "samtools index produced no BAI for ${meta.id}" >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    printf '"%s":\n    samtools: 1.20\n' '${task.process}' > versions.yml
    """
}
