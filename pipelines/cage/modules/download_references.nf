process DOWNLOAD_REFERENCES {
    tag "${meta.species}:${meta.assembly}"
    container 'quay.io/biocontainers/python:3.12.12'
    label 'cage_python'
    publishDir { "${params.outdir}/references/${meta.assembly}" }, mode: 'copy'

    input:
    tuple val(meta), val(fasta_source), val(gtf_source), val(chrom_sizes_source), val(fasta_checksum), val(gtf_checksum), val(chrom_sizes_checksum)
    path helper, stageAs: 'cache_reference.py'

    output:
    tuple val(meta), path('reference.fa'), path('reference.gtf'), path('chrom.sizes'), emit: reference

    script:
    """
    set -euo pipefail
    python cache_reference.py --assembly ${meta.assembly} --cache-dir ${params.reference_cache} \\
        --fasta-source '${fasta_source}' --gtf-source '${gtf_source}' --chrom-sizes-source '${chrom_sizes_source}' \\
        --fasta-checksum '${fasta_checksum}' --gtf-checksum '${gtf_checksum}' --chrom-sizes-checksum '${chrom_sizes_checksum}'
    """

    stub:
    """
    touch reference.fa reference.gtf chrom.sizes
    """
}
