// modules/local/render_stable_id_sql.nf
process RENDER_STABLE_ID_SQL {
    // tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/sql"
    }, mode: 'copy',
    saveAs: { name ->
        name == "${db_name}.stable_id_updates.sql" ||
        name == "${db_name}.stable_id_updates.dry_run.sql" ? name : null
    }

    input:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              path(locus_comparison),
              path(decisions_tsv),
              path(score_evidence_tsv),
              path(decisions_json)

    output:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              path(locus_comparison),
              path(decisions_tsv),
              path(score_evidence_tsv),
              path("${db_name}.stable_id_updates.sql"),
              path("${db_name}.stable_id_updates.dry_run.sql"),
              emit: rendered

    script:
		def translation_flag = params.include_translations ? '' : '--no-translations'
		def replace_flag = params.replace_events_for_session ? '--replace-events-for-session' : ''
		
		"""
		python3 ${projectDir}/bin/render_stable_id_sql.py \
		    --decisions-tsv ${decisions_tsv} \
		    --output-sql ${db_name}.stable_id_updates.sql \
		    --db-name ${db_name} \
		    --mapping-session-id ${mapping_session_id} \
		    --batch-size ${params.batch_size} \
		    ${translation_flag} \
		    ${replace_flag}
		
		python3 ${projectDir}/bin/render_stable_id_sql.py \
		    --decisions-tsv ${decisions_tsv} \
		    --output-sql ${db_name}.stable_id_updates.dry_run.sql \
		    --db-name ${db_name} \
		    --mapping-session-id ${mapping_session_id} \
		    --batch-size ${params.batch_size} \
		    --dry-run \
		    ${translation_flag} \
		    ${replace_flag}
		"""
}
