process REASSIGN_STABLE_IDS {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/reassignment"
    }, mode: 'copy',
    saveAs: { name ->
        name.startsWith("${db_name}.stable_id_reassignment.")
            ? name
            : null
    }

    input:
    tuple val(db_name),
          path(target_gff),
          val(gene_range),
          val(transcript_range),
          val(translation_range)

    output:
    tuple val(db_name),
          path(target_gff, includeInputs: true),
          path("${db_name}.stable_id_reassignment.sql"),
          path("${db_name}.stable_id_reassignment.dry_run.sql"),
          path("${db_name}.stable_id_reassignment.json"),
          path("${db_name}.stable_id_reassignment.id_map.tsv"),
          emit: reassignment

    script:
    """
    python3 ${projectDir}/bin/generate_stable_id_reassignment.py \
        --db-name ${db_name} \
        --gene-range '${gene_range}' \
        --transcript-range '${transcript_range}' \
        --translation-range '${translation_range}' \
        --output-sql ${db_name}.stable_id_reassignment.sql \
        --dry-run-sql ${db_name}.stable_id_reassignment.dry_run.sql \
        --output-json ${db_name}.stable_id_reassignment.json \
        --output-id-map ${db_name}.stable_id_reassignment.id_map.tsv \
        --batch-size ${params.batch_size}
    """
}
