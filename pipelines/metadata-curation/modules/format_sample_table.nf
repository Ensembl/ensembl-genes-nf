process FORMAT_SAMPLE_TABLE {
    tag "$meta.id"
    label 'process_low'

    // No conda/container - use base environment

    input:
    tuple val(meta), path(processed_metadata)
    path biogroups
    val output_format
    val include_raw

    output:
    path "sample_table.*"    , emit: table
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def biogroups_arg = biogroups.name != 'NO_FILE' ? "--biogroups ${biogroups}" : ''
    def raw_arg = include_raw ? "--include-raw" : ''
    def ext = output_format == 'tsv' ? 'tsv' : 'csv'
    """
    python ${projectDir}/bin/table_formatter.py \\
        ${processed_metadata} \\
        -o sample_table.${ext} \\
        ${biogroups_arg} \\
        --format ${output_format} \\
        --level sample \\
        ${raw_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def ext = output_format == 'tsv' ? 'tsv' : 'csv'
    """
    echo "sample_id,study_id,organism" > sample_table.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}