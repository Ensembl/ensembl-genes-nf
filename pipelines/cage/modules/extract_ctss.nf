process EXTRACT_CTSS {
    tag "${meta.sample_id}"
    container 'quay.io/biocontainers/samtools:1.25--h9dcdb79_0'
    label 'cage_samtools'
    publishDir { "${params.outdir}/ctss/${meta.species}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(alignment_log), path(excluded_contigs)
    path helper

    output:
    tuple val(meta), path('*.ctss.bed'), path('*.ctss.bedgraph'), emit: ctss

    script:
    """
    bash ${helper} --bam ${bam} --sample ${meta.sample_id} --mode ${meta.ctss_mode} --exclude-contigs ${excluded_contigs} --bed ${meta.sample_id}.ctss.bed --bedgraph ${meta.sample_id}.ctss.bedgraph
    """

    stub:
    """
    printf 'chr1\\t0\\t1\\t${meta.sample_id}\\t1\\t+\\n' > ${meta.sample_id}.ctss.bed
    printf 'chr1\\t0\\t1\\n' > ${meta.sample_id}.ctss.bedgraph
    """
}
