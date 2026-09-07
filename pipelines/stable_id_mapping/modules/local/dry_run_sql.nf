// modules/local/dry_run_sql.nf
process DRY_RUN_SQL {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/dry_run_sql"
    }, mode: 'copy',
    saveAs: { name -> name.startsWith("${db_name}.") ? name : null }

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
        path("${db_name}.dry_run_sql.out"),
              emit: dry_run_sql

    script:
    """
    gb1-w < ${dry_run_sql} \
        > ${db_name}.dry_run_sql.out
    """
}

