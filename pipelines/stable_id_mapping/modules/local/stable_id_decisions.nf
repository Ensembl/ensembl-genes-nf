// modules/local/stable_id_decisions.nf
process STABLE_ID_DECISIONS {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/decisions"
    }, mode: 'copy',
    saveAs: { name ->
        name == "${db_name}.stable_id_decisions.tsv" ||
        name == "${db_name}.score_evidence.tsv" ||
        name == "${db_name}.stable_id_decisions.json" ? name : null
    }

    input:
        tuple val(db_name),
              path(ref_fasta),
              path(ref_gff),
              path(target_fasta),
              path(target_gff),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range),
              path(projected_gff),
              path(missing_report),
              path(transcript_pairs),
              path(gene_pairs),
              path(locus_comparison),
              path(matching_json)
        path rules_config

    output:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              path(locus_comparison),
              path("${db_name}.stable_id_decisions.tsv"),
              path("${db_name}.score_evidence.tsv"),
              path("${db_name}.stable_id_decisions.json"),
              emit: decisions

    script:
    def translation_flag = params.include_translations ? '' : '--no-translations'
    """
    python3 ${projectDir}/bin/stable_id_decisions.py \
        --db-name ${db_name} \
        --ref-gff ${ref_gff} \
        --target-gff ${target_gff} \
		--ref-fasta ${ref_fasta} \
        --target-fasta ${target_fasta} \
        --mapped-gff ${projected_gff} \
        --missing-report ${missing_report} \
        --transcript-pairs ${transcript_pairs} \
        --gene-pairs ${gene_pairs} \
        --mapping-session-id ${mapping_session_id} \
        --gene-range ${gene_range} \
        --transcript-range ${transcript_range} \
        --translation-range ${translation_range} \
        --decisions-tsv ${db_name}.stable_id_decisions.tsv \
        --score-evidence-tsv ${db_name}.score_evidence.tsv \
        --rules-config ${rules_config} \
        --output-json ${db_name}.stable_id_decisions.json \
        ${translation_flag}
    """
}
