process VERIFY_SOURCE_CHECKSUM {
    tag { "${meta.gca_accession}:${meta.source_entity_id}:${meta.track_type}" }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), path(track)

    output:
    tuple val(meta), path(track), emit: checked
    path "${meta.safe_source_id}.${meta.track_type}.source_validation.tsv", emit: validation
    path "${meta.safe_source_id}.${meta.track_type}.checksum.versions.yml", emit: versions

    script:
    def expected = meta.expected_sha256 ?: ''
    """
    test -s ${track} || { echo "Source track is empty: ${track}" >&2; exit 1; }
    actual=\$(sha256sum ${track} | awk '{print \$1}')
    if [[ -n '${expected}' && '${expected}' != 'stub' && '${expected}' != "\$actual" ]]; then
        echo "Source checksum mismatch for ${meta.source_entity_id}:${meta.track_type}: expected ${expected}, observed \$actual" >&2
        exit 1
    fi
    printf 'source_entity_id\\ttrack_type\\tstatus\\tsha256\\n${meta.source_entity_id}\\t${meta.track_type}\\tcomplete\\t%s\\n' "\$actual" > ${meta.safe_source_id}.${meta.track_type}.source_validation.tsv
    printf '"%s":\\n    sha256sum: repository-shell\\n' '${task.process}' > ${meta.safe_source_id}.${meta.track_type}.checksum.versions.yml
    """

    stub:
    """
    printf 'source_entity_id\\ttrack_type\\tstatus\\tsha256\\n${meta.source_entity_id}\\t${meta.track_type}\\tcomplete\\tstub\\n' > ${meta.safe_source_id}.${meta.track_type}.source_validation.tsv
    printf '"%s":\\n    sha256sum: stub\\n' '${task.process}' > ${meta.safe_source_id}.${meta.track_type}.checksum.versions.yml
    """
}
