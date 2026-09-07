process RESOLVE_SPECIES_INPUTS {
    tag "$db_name ($requested_mode) ($params.assembly_metadata_db)"

	// debug true

	publishDir {
        "${params.output_dir}/${db_name}/inputs"
        }, mode: 'copy',
        saveAs: { name -> name.startsWith("${db_name}.") ? name : null }

    input:
    tuple val(db_name),
          val(requested_mode),
          val(target_fasta),
          val(target_gff),
          val(mapping_session_id),
          val(ref_fasta),
          val(ref_gff)

    output:
    path "${db_name}.species_inputs.json",
         emit: inputs_json

    script:
    """
    python3 ${projectDir}/bin/resolve_species_inputs.py \
        --db-name ${db_name} \
		--assembly-metadata-db '${params.assembly_metadata_db}' \
        --mode ${requested_mode} \
        --target-fasta '${target_fasta}' \
        --target-gff '${target_gff}' \
        --mapping-session-id '${mapping_session_id}' \
        --ref-fasta '${ref_fasta}' \
        --ref-gff '${ref_gff}' \
        --output-json ${db_name}.species_inputs.json
    """
}
