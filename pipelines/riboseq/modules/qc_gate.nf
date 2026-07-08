process QC_GATE {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10 conda-forge::pyyaml"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    publishDir "${params.outdir}/qc_gate", mode: 'copy', pattern: "*.{offsets.pass.tsv,offsets.selected.tsv,pass_lengths.tsv,qc.json,qc_eval.tsv,qc_rule_set.tsv,gate_selection.tsv,translon.selected.txt}"

    input:
    val run_id
    tuple val(meta), path(best_offset), path(metrics)
    path rules

    output:
    tuple val(meta), path("*.offsets.pass.tsv"), emit: filtered_offsets
    tuple val(meta), path("*.offsets.selected.tsv"), optional: true, emit: selected_offsets
    tuple val(meta), path("*.translon.selected.txt"), optional: true, emit: translon_selected
    tuple val(meta), path("*.pass_lengths.tsv"), emit: pass_lengths
    tuple val(meta), path("*.qc.json"), emit: qc_json
    tuple val(meta), path("*.qc_eval.tsv"), emit: qc_eval
    tuple val(meta), path("*.qc_rule_set.tsv"), emit: qc_rule_set
    tuple val(meta), path("*.gate_selection.tsv"), emit: gate_selection

    script:
    def apply_for = params.apply_gate_for ?: 'translon,trackhub'
    def rule_set_name = params.qc_rule_set_name ?: 'default'
    def prefix = meta.id
    """
    qc_gate.py \
      --run-id ${run_id} \
      --sample-id ${meta.id} \
      --metrics-tsv ${metrics} \
      --best-offset ${best_offset} \
      --rules-yaml ${rules} \
      --rule-set-name ${rule_set_name} \
      --apply-for ${apply_for} \
      --out-prefix ${prefix}
    """
}
