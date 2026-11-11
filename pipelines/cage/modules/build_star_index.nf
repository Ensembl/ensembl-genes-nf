process BUILD_STAR_INDEX {
    tag "${meta.assembly}"
    container 'quay.io/biocontainers/star:2.7.11b--h5814d7d_9'
    label 'cage_star'
    publishDir { "${params.outdir}/references/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(gtf), path(chrom_sizes)

    output:
    tuple val(meta), path('star_index'), emit: index

    script:
    """
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf} --runThreadN ${task.cpus}
    """

    stub:
    """
    mkdir star_index
    touch star_index/Genome
    """
}
