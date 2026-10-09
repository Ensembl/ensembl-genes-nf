process REWRITE_TARGET_GFF3 {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/gff3"
    }, mode: 'copy',
    saveAs: { name ->
        name == "${db_name}.stable_ids.gff3" ? name : null
    }

    input:
    tuple val(db_name),
          path(target_gff),
          path(id_map)

    output:
    tuple val(db_name),
          path("${db_name}.stable_ids.gff3"),
          emit: rewritten_gff

    script:
    """
    python3 ${projectDir}/bin/rewrite_target_gff3.py \
        --input-gff ${target_gff} \
        --id-map ${id_map} \
        --output-gff ${db_name}.stable_ids.gff3
    """
}