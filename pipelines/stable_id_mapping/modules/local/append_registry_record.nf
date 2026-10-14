// modules/local/append_registry_record.nf
process APPEND_REGISTRY_RECORD {
	tag "$db_name"

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
              path(executable_sql),
              path(dry_run_sql)

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
    		  emit: appended

	script:
	"""
	append_gene_coverage.py \\
        --decisions-tsv ${decisions_tsv} \\
        --sql-in ${db_name}.incomplete.sql \\
        --sql-out ${db_name}.stable_id_updates.sql \\
        --dry-run-sql-in ${db_name}.incomplete.dry_run.sql \\
        --dry-run-sql-out ${db_name}.stable_id_updates.dry_run.sql \\
        --assembly-metadata-db ${params.assembly_metadata_db} \\
		--mapping-session-id ${mapping_session_id}
	"""
}
