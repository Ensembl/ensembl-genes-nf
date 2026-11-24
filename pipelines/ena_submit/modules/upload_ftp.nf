process ENA_FTP_UPLOAD {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'community.wave.seqera.io/library/curl:8.18.0--78f80c4b644630b0'
    secret 'ENA_WEBIN_PASSWORD'

    input:
    tuple val(meta), val(row), val(file_meta), path(file), path(md5)
    val remote_dir
    val ftp_host
    val webin_user

    output:
    tuple val(meta), val(row), val(file_meta), path(file, includeInputs: true), path(md5, includeInputs: true), emit: uploaded
    path 'versions.yml', emit: versions

    script:
    def basename    = file.getName()
    def remoteName  = file_meta.remote_name ?: basename
    def remoteDir   = remote_dir ? remote_dir.replaceAll('/+$','') : ''
    def remote_path = remoteDir ? "${remoteDir}/${remoteName}" : remoteName
    """
    set -euo pipefail
    command -v curl >/dev/null 2>&1 || { echo 'curl is required' >&2; exit 127; }
    curl --fail --show-error --ftp-pasv --ftp-create-dirs \
        --user "${webin_user}:\$ENA_WEBIN_PASSWORD" \
        --upload-file "${file}" \
        "ftp://${ftp_host}/${remote_path}"
    curl --fail --show-error --ftp-pasv --ftp-create-dirs \
        --user "${webin_user}:\$ENA_WEBIN_PASSWORD" \
        --upload-file "${md5}" \
        "ftp://${ftp_host}/${remote_path}.md5"
    printf 'ENA_FTP_UPLOAD:\n  curl: "%s"\n' "\$(curl --version | awk 'NR==1 {print \$2}')" > versions.yml
    """

    stub:
    """
    touch ${md5}
    printf 'ENA_FTP_UPLOAD:\n  curl: "stub"\n' > versions.yml
    """
}
