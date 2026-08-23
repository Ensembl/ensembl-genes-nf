/* Publish the final indexed STAR alignments as a stable downstream contract.
 * The alignment and indexing tasks can keep using work directories internally,
 * but downstream pipelines need both BAM types and their indexes. */

process PUBLISH_STAR_ALIGNMENTS {
    tag "${meta.id}:${bam_type}"
    label 'process_low'

    publishDir "${params.outdir}/star_align", mode: 'copy', saveAs: { filename ->
        "${bam_type}/${meta.id}/${filename}"
    }

    input:
    tuple val(meta), path(bam), path(bai), val(bam_type)

    output:
    tuple val(meta), path('*.bam'), path('*.bai'), emit: published

    script:
    def suffix = bam_type == 'genome'
        ? 'Aligned.sortedByCoord.out'
        : 'Aligned.toTranscriptome.out'
    """
    if [ "${bam}" != "${meta.id}.${suffix}.bam" ]; then
        cp -L ${bam} ${meta.id}.${suffix}.bam
    fi
    if [ "${bai}" != "${meta.id}.${suffix}.bam.bai" ]; then
        cp -L ${bai} ${meta.id}.${suffix}.bam.bai
    fi
    """

    stub:
    def suffix = bam_type == 'genome'
        ? 'Aligned.sortedByCoord.out'
        : 'Aligned.toTranscriptome.out'
    """
    touch ${meta.id}.${suffix}.bam
    touch ${meta.id}.${suffix}.bam.bai
    """
}
