process WRITE_TRACK_REPORTS {
    tag 'track-report'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path result_files
    path report_script
    path version_files

    output:
    path 'track_manifest.tsv', emit: track_manifest
    path 'entity_status.tsv', emit: entity_status
    path 'exceptions.tsv', emit: exceptions
    path 'versions.yml', emit: versions

    script:
    """
    python ${report_script} --result ${result_files.join(' ')} --output track_manifest.tsv
    cat ${version_files.join(' ')} > versions.yml
    """

    stub:
    """
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > track_manifest.tsv
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\tstatus\n' > entity_status.tsv
    printf 'gca_accession\tassembly_release\tentity_id\tentity_type\ttrack_type\ttrack_path\tsource_path\tstatus\tsha256\ttool_versions\tnormalization_parameters\n' > exceptions.tsv
    printf '"WRITE_TRACK_REPORTS":\n    python: stub\n' > versions.yml
    """
}
