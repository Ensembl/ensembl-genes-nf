// modules/local/dry_run_sql.nf
process DRY_RUN_SQL {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/dry_run_sql"
    }, mode: 'copy',
    saveAs: { name -> name.startsWith("${db_name}.") ? name : null }

    input:
        tuple val(db_name),
              path(dry_run_sql)

    output:
        path("${db_name}.dry_run_sql.out"),
              emit: dry_run_sql

    script:
    """
    mysql -h '${params.gb_host}' -P ${params.gb_port} -u ensadmin -p'${params.ensadmin_password}' \
		< ${dry_run_sql} \
        > ${db_name}.dry_run_sql.out
    """
}

