process ALIGNMENT_QC {
    tag "${meta.id}:alignment-qc"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('alignment_stats.tsv'), emit: stats
    path 'versions.yml', emit: versions

    script:
    """
    samtools quickcheck -v "${bam}" || { echo "Alignment quickcheck failed for ${meta.id}" >&2; exit 1; }
    samtools flagstat "${bam}" > alignment_stats.tsv
    printf 'run\tclassification\tpreset\tsecondary\tbam\n${meta.id}\t${meta.classification}\t${meta.minimap2_preset ?: params.minimap2_preset}\t${params.secondary_mode}\t${bam}\n' >> alignment_stats.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    printf 'run\tclassification\tpreset\tsecondary\tbam\n${meta.id}\t${meta.classification}\t${meta.minimap2_preset ?: params.minimap2_preset}\tstub\t${bam}\n' > alignment_stats.tsv
    printf '"%s":\n    samtools: 1.20\n' '${task.process}' > versions.yml
    """
}
