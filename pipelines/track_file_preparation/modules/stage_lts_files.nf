process STAGE_LTS_FILES {
    tag { "${meta.gca_accession}:${meta.id}" }
    label 'process_light'
    errorStrategy 'retry'
    maxRetries 2
    input:
    tuple val(meta), val(source_fields), val(stage_kinds), val(stage_names), path(sources, stageAs: 'source??')

    output:
    tuple val(meta), val(source_fields), val(stage_kinds), path('stage_*'), emit: staged
    path "${meta.safe_id}.stage.versions.yml", emit: versions

    script:
    def copy_commands = sources.withIndex().collect { source, index ->
        "cp -p -- '${source}' '${stage_names[index]}' && cp -p -- '${source}' \"${params.lts_stage_dir}/${stage_names[index]}\""
    }.join('\n')
    """
    set -euo pipefail
    mkdir -p "${params.lts_stage_dir}"
    ${copy_commands}
    for staged_file in stage_*; do
        test -s "\$staged_file" || { echo "Staged input is empty: \$staged_file" >&2; exit 1; }
    done
    printf '"%s":\\n    coreutils: %s\\n' '${task.process}' "\$(cp --version | head -n1)" > ${meta.safe_id}.stage.versions.yml
    """

    stub:
    def copy_commands = sources.withIndex().collect { source, index ->
        "cp -p -- '${source}' '${stage_names[index]}' && cp -p -- '${source}' \"${params.lts_stage_dir}/${stage_names[index]}\""
    }.join('\n')
    """
    set -euo pipefail
    mkdir -p "${params.lts_stage_dir}"
    ${copy_commands}
    printf '"%s":\\n    coreutils: stub\\n' '${task.process}' > ${meta.safe_id}.stage.versions.yml
    """
}
