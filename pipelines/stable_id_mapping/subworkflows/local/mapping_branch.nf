include { STAGE_SPECIES_INPUTS } from '../../modules/local/stage_species_inputs.nf'
include { LIFTON_PROJECTION } from '../../modules/local/lifton_projection.nf'
include { STRUCTURAL_MATCHING } from '../../modules/local/structural_matching.nf'
include { STABLE_ID_DECISIONS } from '../../modules/local/stable_id_decisions.nf'
include { RENDER_STABLE_ID_SQL } from '../../modules/local/render_stable_id_sql.nf'
include { AUDIT_RUN } from '../../modules/local/audit_run.nf'
include { DRY_RUN_SQL } from '../../modules/local/dry_run_sql.nf'


workflow MAPPING_BRANCH {
    take:
    mapping_ch

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

    RENDER_STABLE_ID_SQL(
        STABLE_ID_DECISIONS.out.decisions
    )

    AUDIT_RUN(
        RENDER_STABLE_ID_SQL.out.rendered
    )

    DRY_RUN_SQL(
	    RENDER_STABLE_ID_SQL.out.rendered.map { db_name, ref_gff, target_gff, mapping_session_id, locus_comparison, decisions_tsv, score_evidence_tsv, executable_sql, dry_run_sql ->
	        tuple(db_name, dry_run_sql)
	    }
	)

    emit:
    audit = AUDIT_RUN.out.audit
    dry_run_sql = DRY_RUN_SQL.out.dry_run_sql
}
