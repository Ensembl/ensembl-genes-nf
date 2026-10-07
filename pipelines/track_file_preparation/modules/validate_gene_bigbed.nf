process VALIDATE_GENE_BIGBED {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(bigbed), path(gtf)

    output:
    tuple val(meta), path(bigbed), path(gtf), emit: validated
    path "${meta.safe_id}.gene_bigbed_validation.tsv", emit: validation

    script:
    """
    test -s ${bigbed} || { echo 'Gene-model BigBed validation failed' >&2; exit 1; }
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gene_bigbed_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gene_bigbed_validation.tsv
    """
}
