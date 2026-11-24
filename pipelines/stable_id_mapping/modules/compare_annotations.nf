process COMPARE_ANNOTATIONS {
    tag "$comparison_name"

    publishDir {
        "${params.output_dir}/${comparison_name}/comparison"
    }, mode: 'copy',
    saveAs: { name ->
        name == "${comparison_name}.relationships.tsv" ? name : null
    }

    input:
        tuple val(comparison_name),
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
        tuple val(comparison_name),
              path("${comparison_name}.relationships.tsv"),
              emit: relationships

    script:
    def translation_flag = params.include_translations ? '' : '--no-translations'
    """
    python3 ${projectDir}/bin/compare_annotations.py \
        --ref-gff ${ref_gff} \
        --target-gff ${target_gff} \
        --mapped-gff ${projected_gff} \
        --missing-report ${missing_report} \
        --transcript-pairs ${transcript_pairs} \
        --gene-pairs ${gene_pairs} \
        --rules-config ${rules_config} \
        --output-tsv ${comparison_name}.relationships.tsv \
        ${translation_flag}
    """
}
