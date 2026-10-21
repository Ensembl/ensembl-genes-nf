process GTF_TO_GENEPRED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred:482--h0b57e2e_0'

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
