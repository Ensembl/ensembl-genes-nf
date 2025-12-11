process CREATE_SAMPLESHEET {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(bed_files)
    
    output:
    tuple val(meta), path("${meta.id}_samplesheet.tsv"), path(bed_files)
    
    script:
    """
    # Create TSV samplesheet with header
    echo -e "tool\\tbed_path" > ${meta.id}_samplesheet.tsv
    
    # Add each tool/bed pair as a row
    for bed_file in ${bed_files}; do
        tool=\$(echo "\${bed_file}" | cut -d'_' -f1)
        echo -e "\${tool}\\t\${bed_file}" >> ${meta.id}_samplesheet.tsv
    done
    
    echo "Samplesheet created for ${meta.id}:" >&2
    cat ${meta.id}_samplesheet.tsv >&2
    
    echo "Files in directory:" >&2
    ls -lh *.bed >&2
    """
}