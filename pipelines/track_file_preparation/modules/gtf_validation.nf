process VALIDATE_GTF {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(gtf), path(chrom_sizes)
    path validator

    output:
    tuple val(meta), path(gtf), path(chrom_sizes), emit: validated
    path "${meta.safe_id}.gtf_validation.tsv", emit: validation

    script:
    """
    python3 "${validator}" "${gtf}"
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gtf_validation.tsv
    """

    stub:
    """
    printf 'entity_id\\ttrack_type\\tstatus\\n${meta.id}\\tgene_model\\tcomplete\\n' > ${meta.safe_id}.gtf_validation.tsv
    """
}
