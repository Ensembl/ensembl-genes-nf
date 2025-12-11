process CREATE_SAMPLESHEET {
    input:
    tuple val(meta), val(tool_names), path(bed_files)
    
    output:
    tuple val(meta), path("${meta.id}_samplesheet.tsv"), path(bed_files)
    
    script:
    """
    # Create TSV samplesheet with header
    echo -e "tool\\tbed_path" > ${meta.id}_samplesheet.tsv
    
    # Add each tool/bed pair as a row
    ${tool_names.withIndex().collect { tool, idx ->
        "echo -e '${tool}\\t${bed_files[idx]}' >> ${meta.id}_samplesheet.tsv"
    }.join('\n')}
    
    echo "Samplesheet created for ${meta.id}:"
    cat ${meta.id}_samplesheet.tsv
    """
}