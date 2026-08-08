process QC_GATE {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.10 conda-forge::pyyaml"
    container "community.wave.seqera.io/library/pip_pyyaml_duckdb_pandas:5ede6677f4262ec2"

    input:
    val run_id
    tuple val(meta), path(best_offset), path(metrics)
    path rules

    output:
    tuple val(meta), path("*.offsets.pass.tsv"), emit: filtered_offsets
    tuple val(meta), path("*.offsets.selected.tsv"), optional: true, emit: selected_offsets
    tuple val(meta), path("*.offsets.good.tsv"), emit: good_offsets
    tuple val(meta), path("*.offsets.great.tsv"), emit: great_offsets
    tuple val(meta), path("*.translon.selected.txt"), optional: true, emit: translon_selected
    tuple val(meta), path("*.pass_lengths.tsv"), emit: pass_lengths
    tuple val(meta), path("*.qc.json"), emit: qc_json
    tuple val(meta), path("*.qc_eval.tsv"), emit: qc_eval
    tuple val(meta), path("*.qc_rule_set.tsv"), emit: qc_rule_set
    tuple val(meta), path("*.gate_selection.tsv"), emit: gate_selection

    script:
    def apply_for = params.apply_gate_for ?: 'translon,trackhub'
    def rule_set_name = params.qc_rule_set_name ?: 'default'
    def prefix = "${meta.id}.${task.index}"
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

    stub:
    """
    touch ${meta.id}.${task.index}.offsets.pass.tsv
    touch ${meta.id}.${task.index}.offsets.selected.tsv
    touch ${meta.id}.${task.index}.offsets.good.tsv
    touch ${meta.id}.${task.index}.offsets.great.tsv
    touch ${meta.id}.${task.index}.pass_lengths.tsv
    touch ${meta.id}.${task.index}.qc.json
    touch ${meta.id}.${task.index}.qc_eval.tsv
    touch ${meta.id}.${task.index}.qc_rule_set.tsv
    touch ${meta.id}.${task.index}.gate_selection.tsv
    """
}
