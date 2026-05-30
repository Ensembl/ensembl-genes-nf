process QC_GATE {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10 conda-forge::duckdb conda-forge::pyyaml"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    publishDir "${params.outdir}/qc_gate", mode: 'copy', pattern: "*.{offsets.pass.tsv,pass_lengths.tsv,qc.json}"

    input:
    val run_id
    tuple val(meta), path(best_offset)
    path rules

    output:
    tuple val(meta), path("*.offsets.pass.tsv"), emit: filtered_offsets
    tuple val(meta), path("*.pass_lengths.tsv"), emit: pass_lengths
    tuple val(meta), path("*.qc.json"), emit: qc_json

    script:
    def db = params.metrics_db ?: "${params.outdir}/pipeline_info/qc.duckdb"
    def apply_for = params.apply_gate_for ?: 'translon,trackhub'
    def rule_set_name = params.qc_rule_set_name ?: 'default'
    def prefix = meta.id
    """
    qc_gate.py \
      --db ${db} \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --best-offset ${best_offset} \
      --rules-yaml ${rules} \
      --rule-set-name ${rule_set_name} \
      --apply-for ${apply_for} \
      --out-prefix ${prefix}
    """
}
