process SORT_GENE_MODELS {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bed12), path(chrom_sizes), path(gtf)

    output:
    tuple val(meta), path("${meta.safe_id}.sorted.bed12"), path(chrom_sizes), path(gtf), emit: sorted

    script:
    """
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n ${bed12} > ${meta.safe_id}.sorted.bed12
    """

    stub:
    """
    cp ${bed12} ${meta.safe_id}.sorted.bed12
    """
}
