process MINIMAP2_TO_SORTED_BAM {
    tag "${meta.id}"
    label 'process_high'
    container 'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0'

    input:
    tuple val(meta), path(reads)
    path minimap_index
    path reference

    output:
    tuple val(meta), path('*.sorted.bam'), path('*.sorted.bam.bai'), emit: bam
    tuple val(meta), path('alignment_stats.tsv'), emit: stats
    path 'versions.yml', emit: versions

    script:
    def secondary = task.ext.secondary ?: (params.secondary_mode == 'yes' ? '--secondary=yes' : '--secondary=no')
    def preset = meta.minimap2_preset ?: params.minimap2_preset
    """
    set -euo pipefail
    test "${meta.classification}" != "PACBIO_SUBREAD_FASTQ_ONLY"
    test "${meta.classification}" != "PACBIO_SUBREAD_BAM"
    minimap2 -t ${task.cpus} -ax ${preset} ${secondary} ${minimap_index} ${reads} \\
      | samtools sort -@ ${task.cpus} -m ${params.samtools_sort_memory} -O BAM -o ${meta.id}.sorted.bam -
    samtools quickcheck -v ${meta.id}.sorted.bam
    samtools index ${meta.id}.sorted.bam
    samtools flagstat ${meta.id}.sorted.bam > alignment_stats.tsv
    printf 'run\\tclassification\\tpreset\\tsecondary\\tbam\\n${meta.id}\\t${meta.classification}\\t${preset}\\t${secondary}\\t${meta.id}.sorted.bam\\n' >> alignment_stats.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai
    printf 'run\\tclassification\\tpreset\\tsecondary\\tbam\\n${meta.id}\\t${meta.classification}\\t${meta.minimap2_preset ?: params.minimap2_preset}\\tstub\\t${meta.id}.sorted.bam\\n' > alignment_stats.tsv
    printf '"%s":\\n    minimap2: 2.28\\n    samtools: 1.20\\n' '${task.process}' > versions.yml
    """
}
