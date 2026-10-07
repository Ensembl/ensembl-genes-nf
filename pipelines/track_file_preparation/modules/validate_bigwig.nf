process VALIDATE_BIGWIG {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bigwig), path(bam)

    output:
    tuple val(meta), path(bigwig), path(bam), emit: validated
    path "${meta.safe_id}.bigwig_validation.tsv", emit: validation

    script:
    """
    test -s ${bigwig} || { echo 'BigWig validation failed: empty file' >&2; exit 1; }
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tcoverage\\tcomplete\\n' > ${meta.safe_id}.bigwig_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tcoverage\\tcomplete\\n' > ${meta.safe_id}.bigwig_validation.tsv
    """
}
