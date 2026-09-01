// modules/local/lifton_projection.nf
process LIFTON_PROJECTION {
    // tag "$db_name"

	publishDir "${params.output_dir}", mode: 'copy',
        saveAs: { name ->
            name.startsWith("${db_name}.")
                ? "${db_name}/lifton/${name}"
                : null
        }

    input:
        tuple val(db_name),
              path(ref_fasta, stageAs: 'reference_fasta/*'),
              path(ref_gff, stageAs: 'reference_gff/*'),
              path(target_fasta, stageAs: 'target_fasta/*'),
              path(target_gff, stageAs: 'target_gff/*'),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range)

    output:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              val(gene_range),
              val(transcript_range),
              val(translation_range),
              path("${db_name}.projected_ref_on_target.gff3"),
              path("${db_name}.missing_genes.txt"),
              path("${db_name}.lifton_projection.json"),
              emit: projected


    script:
    """
    lifton_projection.py \
        --ref-fasta ${ref_fasta} \
        --ref-gff ${ref_gff} \
        --target-fasta ${target_fasta} \
        --output-gff ${db_name}.projected_ref_on_target.gff3 \
        --missing-report ${db_name}.missing_genes.txt \
        --output-json ${db_name}.lifton_projection.json \
        --threads ${task.cpus} \
        --executable ${params.lifton_executable} \
        --feature-types ${params.lifton_feature_types}
    """
}
