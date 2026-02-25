process ENA_COMPUTE_MD5 {
    label 'process_light'
    tag { meta.id }

    input:
    tuple val(meta), val(row), path(file)

    output:
    tuple val(meta), val(row), path(file), path('md5.txt'), emit: md5

    shell:
    """
    md5sum ${file} | awk '{print $1}' > md5.txt
    """
}
