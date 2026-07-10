process ENA_COMPUTE_MD5 {
    label 'process_light'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/library/ubuntu:22.04'

    input:
    tuple val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(meta), val(row), val(file_meta), path(file), path('*.md5'), emit: md5

    shell:
    """
    md5sum "${file}" > "!{file_meta.remote_name ?: file.getName()}.md5"
    """
}
