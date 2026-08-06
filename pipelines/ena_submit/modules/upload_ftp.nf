process ENA_FTP_UPLOAD {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/minidocks/lftp@sha256:0eb637d5fd61bfef82d94a9774880bff9b48e3ff35ac79f0bc5219e095592a94'
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
    def mkdir_cmd   = remoteDir ? "mkdir -p ${remoteDir}; " : ''
    """
    set -euo pipefail
    command -v lftp >/dev/null 2>&1 || { echo 'lftp is required' >&2; exit 127; }
    LFTP_PASSWORD="\$ENA_WEBIN_PASSWORD" \
    lftp --user "${webin_user}" --env-password "${ftp_host}" \
        -e '${mkdir_cmd}put ${file} -o ${remote_path}; put ${md5} -o ${remote_path}.md5; bye'
    printf 'ENA_FTP_UPLOAD:\n  lftp: "%s"\n' "\$(lftp --version | awk 'NR==1 {print \$2}')" > versions.yml
    """

    stub:
    """
    touch ${md5}
    printf 'ENA_FTP_UPLOAD:\n  lftp: "stub"\n' > versions.yml
    """
}
