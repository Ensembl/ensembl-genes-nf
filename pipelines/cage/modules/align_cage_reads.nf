process ALIGN_CAGE_READS {
    tag "${meta.sample_id}:${meta.assembly}"
    label 'cage_star'
    publishDir { "${params.outdir}/alignments/${meta.species}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads), path(index), path(excluded_contigs, stageAs: 'excluded_contigs.input.txt')

    output:
    tuple val(meta), path('*.bam'), path('*.Log.final.out'), path('excluded_contigs.txt'), emit: bam

    script:
    def files = reads.collect { read -> read.toString() }.join(' ')
    """
    set -euo pipefail
    cp excluded_contigs.input.txt excluded_contigs.txt
    STAR --genomeDir ${index} --readFilesIn ${files} --readFilesCommand cat --runThreadN ${task.cpus} ${meta.alignment_args ?: ''} --outFileNamePrefix ${meta.sample_id}. --outSAMtype BAM SortedByCoordinate
    """

    stub:
    """
    touch ${meta.sample_id}.bam ${meta.sample_id}.Log.final.out excluded_contigs.txt
    """
}
