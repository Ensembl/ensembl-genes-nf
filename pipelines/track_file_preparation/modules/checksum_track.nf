process CHECKSUM_TRACK {
    tag { "${meta.gca_accession}:${meta.id}:${meta.current_track}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(track), path(source)

    output:
    tuple val(meta), path(track), path(source), path("${meta.safe_id}.${meta.current_track}.sha256"), emit: checksummed

    script:
    """
    sha256sum ${track} > ${meta.safe_id}.${meta.current_track}.sha256
    """

    stub:
    """
    printf 'stub  ${track}\\n' > ${meta.safe_id}.${meta.current_track}.sha256
    """
}
