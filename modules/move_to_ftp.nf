process MOVE_TO_FTP {
    label 'process_low'

    // Submit to datamover queue/partition
    queue 'datamover'

    tag "${meta.id}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    // Switch to genebuild user context before running main script
    beforeScript 'become genebuild'

    input:
    tuple val(meta), path(files)        // Meta map with file info + files to transfer
    val ftp_destination                 // FTP destination path

    output:
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cp ${args} ${files} ${ftp_destination}

    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cp: 9.1
    END_VERSIONS
    """
}