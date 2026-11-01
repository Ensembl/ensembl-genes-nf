process CHECKSUM_TRACK {
    tag { "${meta.gca_accession}:${meta.id}:${meta.current_track}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(track)

    output:
    tuple val(meta), path(track), path("${meta.safe_id}.${meta.current_track}.sha256"), emit: checksummed
    path "${meta.safe_id}.${meta.current_track}.checksum.versions.yml", emit: versions

    script:
    """
    sha256sum ${track} > ${meta.safe_id}.${meta.current_track}.sha256
    printf '"%s":\\n    sha256sum: repository-shell\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.checksum.versions.yml
    """

    stub:
    """
    printf 'stub  ${track}\\n' > ${meta.safe_id}.${meta.current_track}.sha256
    printf '"%s":\\n    sha256sum: stub\\n' '${task.process}' > ${meta.safe_id}.${meta.current_track}.checksum.versions.yml
    """
}
