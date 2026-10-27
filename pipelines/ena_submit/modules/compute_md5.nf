process ENA_COMPUTE_MD5 {
    label 'process_light'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/library/ubuntu@sha256:0199853f6d6b20b0424f3c5694a72a62764f01e6a771b1eb48a4197848986c7e'

    input:
    tuple val(meta), val(row), val(file_meta), path(file)

    output:
    tuple val(meta), val(row), val(file_meta), path(file), path('*.md5'), emit: md5
    path 'versions.yml', emit: versions

    shell:
    """
    md5sum "${file}" > "!{file_meta.remote_name ?: file.getName()}.md5"
    printf 'ENA_COMPUTE_MD5:\n  coreutils: "%s"\n' "\$(md5sum --version | awk 'NR==1 {print \$NF}')" > versions.yml
    """

    stub:
    """
    touch ${file} ${file}.md5
    printf 'ENA_COMPUTE_MD5:\n  coreutils: "stub"\n' > versions.yml
    """
}
