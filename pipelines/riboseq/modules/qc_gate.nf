process QC_GATE {
    tag "${meta.id}"
    label "process_medium"

    container "python:3.10-slim"

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
    python3 -m pip install -q --no-cache-dir duckdb pyyaml >/dev/null 2>&1 || true
    python3 $projectDir/bin/qc_gate.py \
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
