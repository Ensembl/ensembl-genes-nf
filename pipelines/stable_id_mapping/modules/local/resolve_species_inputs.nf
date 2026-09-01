// modules/local/resolve_species_inputs.nf
process RESOLVE_SPECIES_INPUTS {
    // tag "$db_name"

    input:
    	tuple val(db_name),
    	      val(target_fasta),
    	      val(target_gff),
    	      val(mapping_session_id),
    	      val(ref_fasta),
    	      val(ref_gff)

    output:
        path "${db_name}.species_inputs.json", emit: inputs_json, optional: true

    script:
	"""
	python3 ${projectDir}/bin/resolve_species_inputs.py \
	    --db-name ${db_name} \
	    --target-fasta '${target_fasta}' \
	    --target-gff '${target_gff}' \
	    --mapping-session-id '${mapping_session_id}' \
	    --ref-fasta '${ref_fasta}' \
	    --ref-gff '${ref_gff}' \
	    --output-json ${db_name}.species_inputs.json
	"""
}
