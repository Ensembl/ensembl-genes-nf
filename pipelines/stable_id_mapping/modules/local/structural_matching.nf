// modules/local/structural_matching.nf
process STRUCTURAL_MATCHING {
    // tag "$db_name"

    publishDir {
	    "${params.output_dir}/${db_name}/matching"
	}, mode: 'copy',
	saveAs: { name ->
	    def base = name.tokenize('/').last()
	    base.startsWith("${db_name}.lifton.") ||
	    base == "${db_name}.structural_matching.json" ? base : null
	}

    input:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range),
              path(projected_gff),
              path(missing_report),
              path(projection_json)
        path rules_config

    output:
    	tuple val(db_name),
    	      path(ref_gff),
    	      path(target_gff),
    	      val(mapping_session_id),
    	      val(gene_range),
    	      val(transcript_range),
    	      val(translation_range),
    	      path(projected_gff),
    	      path(missing_report),
    	      path("matching/${db_name}.lifton.transcript_pairs.tsv"),
    	      path("matching/${db_name}.lifton.gene_pairs.tsv"),
    	      path("matching/${db_name}.lifton.gene_locus_comparison.tsv"),
    	      path("${db_name}.structural_matching.json"),
    	      emit: matches

    script:
    """
    structural_matching.py \
        --lifton-gff ${projected_gff} \
        --target-gff ${target_gff} \
        --out-prefix ./matching/${db_name}.lifton \
        --rules-config ${rules_config} \
        --output-json ${db_name}.structural_matching.json
    """
}
