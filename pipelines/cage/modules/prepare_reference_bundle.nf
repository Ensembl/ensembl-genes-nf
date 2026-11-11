process PREPARE_REFERENCE_BUNDLE {
    tag "${meta.species}:${meta.assembly}"
    container 'quay.io/biocontainers/python:3.12.12'
    label 'cage_python'
    publishDir { "${params.outdir}/references/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(gtf), path(chrom_sizes)
    path helper

    output:
    tuple val(meta), path('reference.fa'), path('reference.gtf'), path('chrom.sizes'), path('rrna.bed'), path('excluded_contigs.txt'), emit: reference

    script:
    """
    set -euo pipefail
    cp -L ${gtf} reference.gtf
    cp -L ${chrom_sizes} chrom.sizes
    python ${helper} --fasta ${fasta} --gtf ${gtf} --assembly ${meta.assembly} \\
        --output-fasta reference.fa --output-bed rrna.bed --output-contigs excluded_contigs.txt
    test -s reference.fa || { echo 'Prepared reference FASTA is empty' >&2; exit 1; }
    test -s reference.gtf || { echo 'Reference annotation is empty' >&2; exit 1; }
    test -s chrom.sizes || { echo 'Chromosome sizes are empty' >&2; exit 1; }
    """

    stub:
    """
    touch reference.fa reference.gtf chrom.sizes rrna.bed excluded_contigs.txt
    """
}
