process PREPARE_CORE_RENAME {
    label 'process_light'

    tag "${meta.id}:core-rename"

    publishDir { "${params.outdir}/core_rename" }, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(core_metadata)

    output:
    tuple val(meta), path("core_rename_${meta.id}.csv"), emit: csv
    path 'versions.yml', emit: versions

    script:
    def metadata = new groovy.json.JsonSlurper().parseText(core_metadata.text)
    def oldCoreName = meta.db_name
    if (!(oldCoreName ==~ /.*_1$/)) {
        throw new IllegalArgumentException("Core name does not end in _1: ${oldCoreName}")
    }
    def updateMonth = metadata.last_geneset_update
    def newCoreName = oldCoreName.replaceFirst(/_1$/, "_${updateMonth.replace('-', '')}")

    """
    printf '%s\n' 'core_name,server,new_core_name' > core_rename_${meta.id}.csv
    printf '%s\n' '${oldCoreName},${meta.core_server},${newCoreName}' >> core_rename_${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        core_rename: generated
    END_VERSIONS
    """

    stub:
    """
    printf '%s\n' 'core_name,server,new_core_name' > core_rename_${meta.id}.csv
    printf '%s\n' 'stub_core,stub_server,stub_core_202601' >> core_rename_${meta.id}.csv
    touch versions.yml
    """
}
