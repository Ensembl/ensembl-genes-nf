include { STAGE_TARGET_GFF } from '../../modules/local/stage_target_gff.nf'
include { REASSIGN_STABLE_IDS } from '../../modules/local/reassign_stable_ids.nf'
include { DRY_RUN_SQL } from '../../modules/local/dry_run_sql.nf'
include { REWRITE_TARGET_GFF3 } from '../../modules/local/rewrite_target_gff3.nf'


workflow REASSIGNMENT_BRANCH {
    take:
    reassignment_ch

    main:
    STAGE_TARGET_GFF(reassignment_ch)

    REASSIGN_STABLE_IDS(STAGE_TARGET_GFF.out.staged_reassignment)

    REWRITE_TARGET_GFF3(
        REASSIGN_STABLE_IDS.out.reassignment.map { db_name, target_gff, executable_sql, dry_run_sql, summary_json, id_map ->
            tuple(db_name, target_gff, id_map)
        }
    )

	DRY_RUN_SQL(
		REASSIGN_STABLE_IDS.out.reassignment.map { db_name, target_gff, executable_sql, dry_run_sql, summary_json, id_map ->
			tuple(db_name, dry_run_sql)
		}
	)

    emit:
    reassignment = REASSIGN_STABLE_IDS.out.reassignment
	dry_run_sql = DRY_RUN_SQL.out.dry_run_sql
    gff3 = REWRITE_TARGET_GFF3.out.rewritten_gff
}
