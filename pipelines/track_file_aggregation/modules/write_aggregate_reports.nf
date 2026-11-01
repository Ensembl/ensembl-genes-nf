process WRITE_AGGREGATE_REPORTS {
    tag 'aggregate-report'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path result_files
    path version_files
    path plan_file
    path report_script

    output:
    path 'track_manifest.tsv', emit: track_manifest
    path 'aggregate_status.tsv', emit: aggregate_status
    path 'exceptions.tsv', emit: exceptions
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${report_script} --result ${result_files.join(' ')} --plan ${plan_file} --output track_manifest.tsv
    cat ${version_files.join(' ')} > versions.yml
    """

    stub:
    """
    python3 ${report_script} --result ${result_files.join(' ')} --plan ${plan_file} --output track_manifest.tsv
    printf '"%s":\\n    aggregate_reports: stub\\n' '${task.process}' > versions.yml
    """
}
