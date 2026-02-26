process ENA_COMPUTE_MD5 {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/ubuntu:22.04'

    input:
    tuple val(meta), val(row), path(file)

    output:
    tuple val(meta), val(row), path(file), path('md5.txt'), emit: md5

    shell:
    """
    md5sum "${file}" > md5.txt
    """
}
