process ENA_FTP_UPLOAD {
    label 'process_single_long'
    tag { meta.id }
    container 'minidocks/lftp'

    input:
    tuple val(meta), val(row), path(file), path(md5)
    val remote_dir
    val ftp_host
    val webin_user
    val webin_password

    output:
    tuple val(meta), val(row), path(file, includeInputs: true), path(md5, includeInputs: true), emit: uploaded

    script:
    def basename    = file.getName()
    def remoteName  = meta.remote_name ?: basename
    def remoteDir   = remote_dir ? remote_dir.replaceAll('/+$','') : ''
    def remote_path = remoteDir ? "${remoteDir}/${remoteName}" : remoteName
    """
    set -euo pipefail
    command -v lftp >/dev/null 2>&1 || { echo 'lftp is required' >&2; exit 127; }
    lftp -u ${webin_user},${webin_password} ${ftp_host} \
        -e 'put ${file} -o ${remote_path}; put ${md5} -o ${remote_path}.md5; bye'
    """
}