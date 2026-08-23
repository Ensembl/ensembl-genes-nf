/* Publish the final indexed STAR alignments as a stable downstream contract.
 * The alignment and indexing tasks can keep using work directories internally,
 * but downstream pipelines need both BAM types and their indexes. */

process PUBLISH_STAR_ALIGNMENTS {
    tag "${meta.id}:${bam_type}"
    label 'process_low'

    publishDir "${params.outdir}/star_align", mode: 'copy', saveAs: { filename ->
        def suffix = bam_type == 'genome'
            ? 'Aligned.sortedByCoord.out'
            : 'Aligned.toTranscriptome.out'
        def published_name = filename.endsWith('.bai')
            ? "${meta.id}.${suffix}.bam.bai"
            : "${meta.id}.${suffix}.bam"
        "${bam_type}/${meta.id}/${published_name}"
    }

    input:
    tuple val(meta), path(bam), path(bai), val(bam_type)

    output:
    tuple val(meta), path('published.bam'), path('published.bam.bai'), emit: published

    script:
    """
    cp -L ${bam} published.bam
    cp -L ${bai} published.bam.bai
    """

    stub:
    """
    touch published.bam
    touch published.bam.bai
    """
}
