process ENA_FTP_UPLOAD {
    label 'process_single_long'
    tag { meta.id }

    input:
    tuple val(meta), val(row), path(file), path(md5)
    val remote_dir
    val ftp_host

    output:
    tuple val(meta), val(row), path(file), path(md5), val(remote_path), emit: uploaded

    shell:
    def basename   = file.getName()
    def remoteName = meta.remote_name ?: basename
    def remoteDir  = remote_dir ? (remote_dir.endsWith('/') ? remote_dir : remote_dir + '/') : ''
    def remote_path = "${remoteDir}${remoteName}"
    """
    set -euo pipefail
    : "${WEBIN_USER:?WEBIN_USER not set}"
    : "${WEBIN_PASSWORD:?WEBIN_PASSWORD not set}"
    command -v lftp >/dev/null 2>&1 || { echo 'lftp is required' >&2; exit 127; }

    lftp -e "set ssl:verify-certificate no; open -u $WEBIN_USER,$WEBIN_PASSWORD ftp://${ftp_host}; mkdir -p ${remoteDir}; put -O ${remoteDir} ${file} -o ${remoteName}; put -O ${remoteDir} ${md5} -o ${remoteName}.md5; bye"
    """
}
