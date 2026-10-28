process STABLE_ID_MAPPING_RUN {
    tag "$db_name"

    input:
        tuple val(db_name), path(ref_fasta), path(ref_gff),
              path(target_fasta), path(target_gff),
              val(mapping_session_id), val(gene_range),
              val(transcript_range), val(translation_range)

    output:
        path "${db_name}_out", emit: run_output

    script:
    def dry_run_flag = params.dry_run_sql ? "--dry-run-sql" : ""
    def replace_flag = params.replace_events_for_session ? "--replace-events-for-session" : ""
    """
    python3 ${projectDir}/../pipelines/stable_id_mapper/run_stable_id_mapping.py \\
        --ref-fasta ${ref_fasta} \\
        --ref-gff ${ref_gff} \\
        --target-fasta ${target_fasta} \\
        --target-gff ${target_gff} \\
        --db-name ${db_name} \\
        --mapping-session-id ${mapping_session_id} \\
        --gene-range ${gene_range} \\
        --transcript-range ${transcript_range} \\
        --translation-range ${translation_range} \\
        --output-dir ${db_name}_out \\
        --rules-config ${params.rules_config} \\
        --lifton-threads ${params.lifton_threads} \\
        --lifton-executable ${params.lifton_executable} \\
        --batch-size ${params.batch_size} \\
        ${dry_run_flag} \\
        ${replace_flag}
    """
}
