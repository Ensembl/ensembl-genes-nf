process AUDIT_NATIVE_MODELS {
    tag "${backend}:${meta.id}:${scope}:native-audit"
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    tuple val(meta), val(backend), val(scope), path(native_model), val(native_format)

    output:
    tuple val(meta), val(backend), val(scope), path("${backend}_${meta.id}_${scope}_native_model_manifest.tsv"), emit: manifest
    tuple val(meta), val(backend), val(scope), path("${backend}_${meta.id}_${scope}_native_model_stats.tsv"), emit: stats
    path "${backend}_${meta.id}_${scope}_native_model.sha256", emit: checksum
    path 'versions.yml', emit: versions

    script:
    """
    audit_native_models.py ${native_model} ${native_format} '${backend}' '${scope}' \\
        ${backend}_${meta.id}_${scope}_native_model_manifest.tsv \\
        ${backend}_${meta.id}_${scope}_native_model_stats.tsv \\
        ${backend}_${meta.id}_${scope}_native_model.sha256
    printf '"%s":\\n    native_model_audit: python\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'backend\\tscope\\tnative_format\\tmodel_path\\tmodel_sha256\\tstatus\\n${backend}\\t${scope}\\t${native_format}\\tstub\\tstub\\tCOMPLETE\\n' > ${backend}_${meta.id}_${scope}_native_model_manifest.tsv
    printf 'backend\\tscope\\tnative_format\\tmodel_count\\tstatus\\n${backend}\\t${scope}\\t${native_format}\\t1\\tCOMPLETE\\n' > ${backend}_${meta.id}_${scope}_native_model_stats.tsv
    printf 'stub  stub\\n' > ${backend}_${meta.id}_${scope}_native_model.sha256
    printf '"%s":\\n    native_model_audit: stub\\n' '${task.process}' > versions.yml
    """
}
