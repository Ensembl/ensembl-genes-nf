process GENEPRED_TO_BED12 {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(genepred), path(chrom_sizes), path(gtf)

    output:
    tuple val(meta), path("${meta.safe_id}.bed12"), path(chrom_sizes), path(gtf), emit: bed12

    script:
    """
    genePredToBed ${genepred} ${meta.safe_id}.bed12
    """

    stub:
    """
    printf 'chr1\\t0\\t1\\tstub\\t0\\t+\\t0\\t1\\t0\\t1\\t1\\t0\\n' > ${meta.safe_id}.bed12
    """
}
