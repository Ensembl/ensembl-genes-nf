process WRITE_SUPPORTING_EVIDENCE_REPORTS {
    tag 'supporting-evidence-reports'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

    input:
    path plan
    path result_files
    path writer

    output:
    path 'supporting_evidence_manifest.tsv', emit: manifest
    path 'supporting_evidence_status.tsv', emit: status
    path 'supporting_evidence_metadata.json', emit: metadata
    path 'supporting_evidence_reports.versions.yml', emit: versions

    script:
    """
    python3 ${writer} --plan ${plan} --files ${result_files.join(' ')}
    printf '"%s":\n    write_supporting_evidence_reports: python\n' '${task.process}' > supporting_evidence_reports.versions.yml
    """

    stub:
    """
    python3 ${writer} --plan ${plan} --files ${result_files.join(' ')}
    printf '"%s":\n    write_supporting_evidence_reports: stub\n' '${task.process}' > supporting_evidence_reports.versions.yml
    """
}
