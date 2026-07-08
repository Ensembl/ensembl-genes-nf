process IMPORT_QC_DB {
    label "process_medium"
    maxForks 1

    conda "conda-forge::python=3.10 conda-forge::duckdb"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: "qc.duckdb"

    input:
    path metrics
    path artifacts
    path rule_sets
    path qc_evals
    path gate_selections

    output:
    path "qc.duckdb", emit: db

    script:
    def metric_files = metrics instanceof List ? metrics : [metrics]
    def artifact_files = artifacts instanceof List ? artifacts : [artifacts]
    def rule_set_files = rule_sets instanceof List ? rule_sets : [rule_sets]
    def qc_eval_files = qc_evals instanceof List ? qc_evals : [qc_evals]
    def gate_selection_files = gate_selections instanceof List ? gate_selections : [gate_selections]
    def metrics_args = metric_files.findAll { it }.collect { it.toString() }.join(' ')
    def artifacts_args = artifact_files.findAll { it }.collect { it.toString() }.join(' ')
    def rule_set_args = rule_set_files.findAll { it }.collect { it.toString() }.join(' ')
    def qc_eval_args = qc_eval_files.findAll { it }.collect { it.toString() }.join(' ')
    def gate_selection_args = gate_selection_files.findAll { it }.collect { it.toString() }.join(' ')
    """
    import_qc_db.py \\
      --db qc.duckdb \\
      --metrics ${metrics_args} \\
      --artifacts ${artifacts_args} \\
      --rule-sets ${rule_set_args} \\
      --qc-evals ${qc_eval_args} \\
      --gate-selections ${gate_selection_args}
    """
}
