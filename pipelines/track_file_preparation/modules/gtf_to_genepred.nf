process GTF_TO_GENEPRED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(gtf), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.safe_id}.genePred"), path(chrom_sizes), path(gtf), emit: genepred

    script:
    """
    gtfToGenePred -genePredExt ${gtf} ${meta.safe_id}.genePred
    """

    stub:
    """
    printf 'stub\\n' > ${meta.safe_id}.genePred
    """
}
