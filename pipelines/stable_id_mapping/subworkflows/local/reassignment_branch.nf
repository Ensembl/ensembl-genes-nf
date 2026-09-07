include { REASSIGN_STABLE_IDS } from '../../modules/local/reassign_stable_ids.nf'
include { DRY_RUN_SQL } from '../../modules/local/dry_run_sql.nf'


workflow REASSIGNMENT_BRANCH {
    take:
    reassignment_ch

    main:
    REASSIGN_STABLE_IDS(reassignment_ch)

	DRY_RUN_SQL(
		REASSIGN_STABLE_IDS.out.reassignment.map { db_name, reassignment_sql, dry_run_sql, json ->
			tuple(db_name, dry_run_sql)
		}
	)

    emit:
    reassignment = REASSIGN_STABLE_IDS.out.reassignment
	dry_run_sql = DRY_RUN_SQL.out.dry_run_sql
}
