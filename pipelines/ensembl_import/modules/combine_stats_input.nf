process COMBINE_STATS_INPUT {
    label 'process_light'

    publishDir { "${params.outdir}/statistics_input" }, mode: 'copy', overwrite: true

    input:
    path csv_files

    output:
    path 'statistics_input.csv', emit: csv
    path 'versions.yml', emit: versions

    script:
    """
    first=true
    for csv in ${csv_files}; do
        if \$first; then
            cat \"\$csv\" > statistics_input.csv
            first=false
        else
            tail -n +2 \"\$csv\" >> statistics_input.csv
        fi
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        statistics_input: combined
    END_VERSIONS
    """

    stub:
    """
    cp ${csv_files[0]} statistics_input.csv
    touch versions.yml
    """
}
