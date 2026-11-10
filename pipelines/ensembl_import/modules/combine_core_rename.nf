process COMBINE_CORE_RENAME {
    label 'process_light'

    publishDir { "${params.outdir}/core_rename" }, mode: 'copy', overwrite: true

    input:
    path csv_files

    output:
    path 'core_rename.csv', emit: csv
    path 'versions.yml', emit: versions

    script:
    """
    first=true
    for csv in ${csv_files}; do
        if \$first; then
            cat "\$csv" > core_rename.csv
            first=false
        else
            tail -n +2 "\$csv" >> core_rename.csv
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        core_rename: combined
    END_VERSIONS
    """

    stub:
    """
    cp ${csv_files[0]} core_rename.csv
    touch versions.yml
    """
}
