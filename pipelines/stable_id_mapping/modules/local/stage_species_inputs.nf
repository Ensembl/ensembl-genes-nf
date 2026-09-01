// modules/local/stage_species_inputs.nf
process STAGE_SPECIES_INPUTS {
    // tag "$db_name"

    publishDir {
    	"${params.output_dir}/${db_name}/inputs"
	}, mode: 'copy',
	saveAs: { name -> name.startsWith("${db_name}.") ? name : null }

    input:
        tuple val(db_name),
              path(ref_fasta, stageAs: 'ref_input/*'),
          	  path(ref_gff, stageAs: 'ref_input/*'),
          	  path(target_fasta, stageAs: 'target_input/*'),
          	  path(target_gff, stageAs: 'target_input/*'),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range)

    output:
        tuple val(db_name),
              path("${db_name}.ref.fa"),
              path("${db_name}.ref.gff3"),
              path("${db_name}.target.fa"),
              path("${db_name}.target.gff3"),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range),
              emit: staged_inputs

    script:
	"""
	copy_or_decompress() {
	    input="\$1"
	    output="\$2"
	
	    case "\$input" in
	        *.gz)
	            gzip -cd "\$input" > "\$output"
	            ;;
	        *)
	            cp -L "\$input" "\$output"
	            ;;
	    esac
	}
	
	copy_or_decompress ${ref_fasta} ${db_name}.ref.fa
	copy_or_decompress ${ref_gff} ${db_name}.ref.gff3
	copy_or_decompress ${target_fasta} ${db_name}.target.fa
	copy_or_decompress ${target_gff} ${db_name}.target.gff3
	"""
}
