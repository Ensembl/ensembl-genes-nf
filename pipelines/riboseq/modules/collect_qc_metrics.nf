process COLLECT_QC_METRICS {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"
    
    input:
    val run_id
    tuple val(meta), path(star_log), path(getrpf_report), path(getrpf_checks), path(ribometric_json), path(ribometric_csv), path(offsets)

    output:
    tuple val(meta), path("*.qc_metrics.tsv"), emit: metrics
    tuple val(meta), path("*.qc_artifacts.tsv"), emit: artifacts

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}.${task.index}"
    """
    collect_qc_metrics.py \\
      --run-id ${run_id} \\
      --sample-id ${meta.id} \\
      --study-id ${meta.study_id ?: 'unknown'} \\
      --star-log ${star_log} \\
      --getrpf-report ${getrpf_report} \\
      --getrpf-checks ${getrpf_checks} \\
      --ribometric-json ${ribometric_json} \\
      --ribometric-csv ${ribometric_csv} \\
      --offsets ${offsets} \\
      --out-prefix ${prefix}
    """
}
