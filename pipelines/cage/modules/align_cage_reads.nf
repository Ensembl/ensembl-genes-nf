process ALIGN_CAGE_READS {
    tag "${meta.sample_id}:${meta.assembly}"
    container 'quay.io/biocontainers/star:2.7.11b--h5814d7d_9'
    label 'cage_star'
    publishDir { "${params.outdir}/alignments/${meta.species}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), path('*.Log.final.out'), emit: bam

    script:
    def files = reads.collect { read -> read.toString() }.join(' ')
    """
    STAR --genomeDir ${index} --readFilesIn ${files} --readFilesCommand cat --runThreadN ${task.cpus} ${meta.alignment_args ?: ''} --outFileNamePrefix ${meta.sample_id}. --outSAMtype BAM SortedByCoordinate
    """

    stub:
    """
    touch ${meta.sample_id}.bam ${meta.sample_id}.Log.final.out
    """
}
