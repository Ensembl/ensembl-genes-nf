process BUILD_STAR_INDEX {
    tag "${meta.assembly}"
    container 'quay.io/biocontainers/star:2.7.11b--h5814d7d_9'
    label 'cage_star'
    publishDir { "${params.outdir}/references/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(gtf), path(chrom_sizes), path(rrna_bed), path(excluded_contigs)

    output:
    tuple val(meta), path('star_index'), path('excluded_contigs.txt'), emit: index

    script:
    """
    set -euo pipefail
    cp ${excluded_contigs} excluded_contigs.txt
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${fasta} --sjdbGTFfile ${gtf} --runThreadN ${task.cpus}
    """

    stub:
    """
    mkdir star_index
    touch star_index/Genome excluded_contigs.txt
    """
}
