process RENAME_BED {
    tag "${meta.id} - ${tool}"

    label 'process_light'
    
    input:
    tuple val(meta), val(tool), path(bed_file)
    
    output:
    tuple val(meta), path("${tool}_${bed_file.name}")
    
    script:
    """
    cp ${bed_file} ${tool}_${bed_file.name}
    """
}