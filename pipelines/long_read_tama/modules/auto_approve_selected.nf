process AUTO_APPROVE_SELECTED_LONG_READS {
    tag 'automatic-selection'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11'

    input:
    path classification_report
    path approver

    output:
    path 'auto_approved_run_manifest.tsv', emit: manifest
    path 'automatic_selection_audit.tsv', emit: audit
    path 'versions.yml', emit: versions

    script:
    """
    python3 ${approver} ${classification_report} auto_approved_run_manifest.tsv automatic_selection_audit.tsv
    printf '"%s":\\n    python: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    """
    printf 'run_accession\ttissue\tdescription\tclassification\tproposed_action\tselected_artifact_uri\tselected_artifact_md5\tselected_artifact_basename\tminimap2_preset\texpected_header_representation\tclassification_report_sha256\treviewer\treviewed_at\tstatus\treview_decision\n' > auto_approved_run_manifest.tsv
    printf 'SRR000001\tstub\tstub\tONT_FASTQ\tALIGN_ONT\thttps://example.org/SRR000001.fastq.gz\td41d8cd98f00b204e9800998ecf8427e\tSRR000001.fastq.gz\tsplice\tONT\tstub\tautomatic_selector\tstub\tAPPROVED\tAPPROVE\n' >> auto_approved_run_manifest.tsv
    printf 'run_accession\tclassification\tstatus\tdecision\tdetail\nSRR000001\tONT_FASTQ\tAPPROVED\tPROCESS\tstub\n' > automatic_selection_audit.tsv
    printf '"%s":\\n    python: stub\\n' '${task.process}' > versions.yml
    """
}
