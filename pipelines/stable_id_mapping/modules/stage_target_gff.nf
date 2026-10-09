process STAGE_TARGET_GFF {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/inputs"
    }, mode: 'copy',
    saveAs: { name ->
        name == "${db_name}.target.gff3" ? name : null
    }

    input:
    tuple val(db_name),
          path(target_gff, stageAs: 'target_input/*'),
          val(gene_range),
          val(transcript_range),
          val(translation_range)

    output:
    tuple val(db_name),
          path("${db_name}.target.gff3"),
          val(gene_range),
          val(transcript_range),
          val(translation_range),
          emit: staged_reassignment

    script:
    """
    case "${target_gff}" in
        *.gz)
            gzip -cd "${target_gff}" > "${db_name}.target.gff3"
            ;;
        *)
            cp -L "${target_gff}" "${db_name}.target.gff3"
            ;;
    esac
    """
}