#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RESOLVE_SPECIES_INPUTS } from './modules/local/resolve_species_inputs.nf'
include { STAGE_SPECIES_INPUTS } from './modules/local/stage_species_inputs.nf'
include { LIFTON_PROJECTION } from './subworkflows/lifton_projection.nf'
include { STRUCTURAL_MATCHING } from './modules/local/structural_matching.nf'
include { STABLE_ID_DECISIONS } from './modules/local/stable_id_decisions.nf'
include { RENDER_STABLE_ID_SQL } from './modules/local/render_stable_id_sql.nf'
include { AUDIT_RUN } from './modules/local/audit_run.nf'

workflow {
    if (params.samplesheet == null) {
        error "Provide --samplesheet <path.csv>"
    }

    species_ch = Channel
    	.fromPath(params.samplesheet)
    	.splitCsv(header: true)
    	.map { row ->
    	    tuple(
    	        row.db_name,
    	        row.get('target_fasta') ?: '',
    	        row.get('target_gff') ?: '',
    	        row.get('mapping_session_id') ?: '',
    	        row.get('ref_fasta') ?: '',
    	        row.get('ref_gff') ?: ''
    	    )
    	}

    RESOLVE_SPECIES_INPUTS(species_ch)

	resolved_ch = RESOLVE_SPECIES_INPUTS.out.inputs_json
    .map { json_file ->
        def data = new groovy.json.JsonSlurper().parse(json_file)
        tuple(
            data.db_name,
            file(data.ref_fasta),
            file(data.ref_gff),
            file(data.target_fasta),
            file(data.target_gff),
            data.mapping_session_id,
            data.gene_range,
            data.transcript_range,
            data.translation_range,
        )
    }

	STAGE_SPECIES_INPUTS(resolved_ch)

	LIFTON_PROJECTION(STAGE_SPECIES_INPUTS.out.staged_inputs)

	rules_config_ch = Channel.value(file(params.rules_config, checkIfExists: true))

	STRUCTURAL_MATCHING(LIFTON_PROJECTION.out.projected, rules_config_ch)

	STABLE_ID_DECISIONS(STRUCTURAL_MATCHING.out.matches, rules_config_ch)

	RENDER_STABLE_ID_SQL(STABLE_ID_DECISIONS.out.decisions)

	AUDIT_RUN(RENDER_STABLE_ID_SQL.out.rendered)
	
	// CLEANUP(copied files, non relevant temporary files and so)
}
