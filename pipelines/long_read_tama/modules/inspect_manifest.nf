process INSPECT_LONG_READ_MANIFEST {
    tag 'read-inventory'
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/python:3.11--he2b4eab_0'

    input:
    path manifest
    path metadata
    path inspector

    output:
    path 'classification_report', emit: reports

    script:
    """
    mkdir -p classification_report
    python3 ${inspector} inspect ${manifest} ${metadata} classification_report
    """

    stub:
    """
    mkdir -p classification_report
    printf 'run_accession\\tclassification\\tconfidence\\tstatus\\treason_codes\\nSRR000001\\tUNKNOWN\\tINSUFFICIENT\\tQUARANTINED\\tSTUB\\n' > classification_report/run_classification.tsv
    printf 'run_accession\\tsource\\turi\\tbasename\\tmd5\\tformat\\tartifact_role\\tselected\\n' > classification_report/run_artifacts.tsv
    printf 'run_accession\\theader_representation\\trecords_sampled\\tsubread_count\\tccs_count\\tont_uuid_count\\tmalformed_count\\tprobe_failure\\tfirst_tokens\\nSRR000001\\tUNKNOWN\\t0\\t0\\t0\\t0\\t0\\tstub\\tunknown\\n' > classification_report/header_probe.tsv
    printf 'run_accession\\trecords_sampled\\tdistinct_sampled_molecules\\tstatus\\nSRR000001\\t0\\t0\\tPROBE_ONLY\\n' > classification_report/molecule_audit.tsv
    printf 'classification\\tconfidence\\truns\\nUNKNOWN\\tINSUFFICIENT\\t1\\n' > classification_report/classification_summary.tsv
    printf 'reason_code\\truns\\nSTUB\\t1\\n' > classification_report/reason_code_summary.tsv
    """
}
