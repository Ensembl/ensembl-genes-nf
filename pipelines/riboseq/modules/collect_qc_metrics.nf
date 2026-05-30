process COLLECT_QC_METRICS {
    tag "${meta.id}"
    label "process_medium"
    maxForks 1

    conda "conda-forge::python=3.10 conda-forge::duckdb conda-forge::pandas conda-forge::pyyaml"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"
    
    publishDir "${params.outdir}/qc_gate", mode: 'copy', pattern: "*.{offsets.pass.tsv,pass_lengths.tsv,qc.json}"

    input:
    val run_id
    tuple val(meta), path(star_log), path(getrpf_report), path(getrpf_checks), path(ribometric_json), path(ribometric_csv), path(offsets)
    path rules

    output:
    tuple val(meta), path("*.offsets.pass.tsv"), emit: filtered_offsets
    tuple val(meta), path("*.pass_lengths.tsv"), emit: pass_lengths
    tuple val(meta), path("*.qc.json"), emit: qc_json

    when:
    task.ext.when == null || task.ext.when

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    def apply_for = params.apply_gate_for ?: 'translon,trackhub'
    def rule_set_name = params.qc_rule_set_name ?: 'default'
    def prefix = meta.id
    """
    collect_star_log.py \\
      --db ${db} \\
      --run-id ${run_id} \\
      --sample-id ${meta.id} \\
      --study-id ${meta.study_id ?: 'unknown'} \\
      --log ${star_log}

    if [ -s ${getrpf_report} ] && [ -s ${getrpf_checks} ]; then
      collect_getrpf_clean.py \\
        --db ${db} \\
        --run-id ${run_id} \\
        --sample-id ${meta.id} \\
        --report ${getrpf_report} \\
        --checks ${getrpf_checks}
    fi

    collect_ribometric.py \\
      --db ${db} \\
      --run-id ${run_id} \\
      --sample-id ${meta.id} \\
      --json ${ribometric_json} \\
      --csv ${ribometric_csv} \\
      --offsets ${offsets}

    qc_gate.py \\
      --db ${db} \\
      --run-id ${run_id} \\
      --sample-id ${meta.id} \\
      --best-offset ${offsets} \\
      --rules-yaml ${rules} \\
      --rule-set-name ${rule_set_name} \\
      --apply-for ${apply_for} \\
      --out-prefix ${prefix}
    """
}
