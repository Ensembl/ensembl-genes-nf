process VALIDATE_SPLICE_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bigbed), path(bam)

    output:
    tuple val(meta), path(bigbed), path(bam), emit: validated
    path "${meta.safe_id}.splice_bigbed_validation.tsv", emit: validation

    script:
    """
    test -s ${bigbed} || { echo 'Splice-junction BigBed validation failed' >&2; exit 1; }
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.splice_bigbed_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tsplice_junction\\tcomplete\\n' > ${meta.safe_id}.splice_bigbed_validation.tsv
    """
}
