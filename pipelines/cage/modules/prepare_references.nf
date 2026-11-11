process PREPARE_REFERENCES {
    tag "${meta.assembly}"
    container 'quay.io/biocontainers/python:3.12.12'
    label 'cage_python'
    publishDir { "${params.outdir}/references/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(gtf), path(chrom_sizes)

    output:
    tuple val(meta), path('reference.fa'), path('reference.gtf'), path('chrom.sizes'), emit: reference

    script:
    """
    cp -L ${fasta} reference.fa
    cp -L ${gtf} reference.gtf
    cp -L ${chrom_sizes} chrom.sizes
    test -s reference.fa || { echo 'Reference FASTA is empty' >&2; exit 1; }
    test -s reference.gtf || { echo 'Reference annotation is empty' >&2; exit 1; }
    test -s chrom.sizes || { echo 'Chromosome sizes are empty' >&2; exit 1; }
    """

    stub:
    """
    touch reference.fa reference.gtf chrom.sizes
    """
}
