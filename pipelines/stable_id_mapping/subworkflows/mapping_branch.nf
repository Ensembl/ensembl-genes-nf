include { STAGE_SPECIES_INPUTS } from '../modules/stage_species_inputs.nf'
include { LIFTON_PROJECTION } from '../modules/lifton_projection.nf'
include { STRUCTURAL_MATCHING } from '../modules/structural_matching.nf'
include { STABLE_ID_DECISIONS } from '../modules/stable_id_decisions.nf'
include { RENDER_STABLE_ID_SQL } from '../modules/render_stable_id_sql.nf'
include { APPEND_REGISTRY_RECORD } from '../modules/append_registry_record.nf'
include { AUDIT_RUN } from '../modules/audit_run.nf'
include { DRY_RUN_SQL } from '../modules/dry_run_sql.nf'
include { REWRITE_TARGET_GFF3 } from '../modules/rewrite_target_gff3.nf'


workflow MAPPING_BRANCH {
	take:
	mapping_ch
	mapping_metadata_ch

	main:
	rules_config_ch = Channel.value(
		file(params.rules_config, checkIfExists: true)
	)

	STAGE_SPECIES_INPUTS(mapping_ch)

	LIFTON_PROJECTION(
		STAGE_SPECIES_INPUTS.out.staged_inputs
	)

	STRUCTURAL_MATCHING(
		LIFTON_PROJECTION.out.projected,
		rules_config_ch
	)

	STABLE_ID_DECISIONS(
		STRUCTURAL_MATCHING.out.matches,
		rules_config_ch
	)

    REWRITE_TARGET_GFF3(
		STABLE_ID_DECISIONS.out.decisions.map { db_name, ref_gff, target_gff, mapping_session_id, locus_comparison, decisions_tsv, score_evidence_tsv, decisions_json ->
            tuple(db_name, target_gff, decisions_tsv)
		}
	)

	render_inputs_ch = STABLE_ID_DECISIONS.out.decisions
        .join(mapping_metadata_ch, by: 0)

	RENDER_STABLE_ID_SQL(
		render_inputs_ch
	)

	APPEND_REGISTRY_RECORD(
		RENDER_STABLE_ID_SQL.out.rendered
	)

	AUDIT_RUN(
		RENDER_STABLE_ID_SQL.out.rendered
	)

	DRY_RUN_SQL(
		APPEND_REGISTRY_RECORD.out.appended.map { db_name, ref_gff, target_gff, mapping_session_id, locus_comparison, decisions_tsv, score_evidence_tsv, executable_sql, dry_run_sql ->
			tuple(db_name, dry_run_sql)
		}
	)

	emit:
	audit = AUDIT_RUN.out.audit
	dry_run_sql = DRY_RUN_SQL.out.dry_run_sql
	gff3 = REWRITE_TARGET_GFF3.out.rewritten_gff
}
