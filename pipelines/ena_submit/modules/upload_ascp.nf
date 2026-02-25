process ENA_ASCP_UPLOAD {
    label 'process_single_long'
    tag { meta.id }

    input:
    tuple val(meta), val(row), path(file), path(md5)
    val remote_dir
    val ascp_limit
    val ascp_host

    output:
    tuple val(meta), val(row), path(file), path(md5), val(remote_path), emit: uploaded

    when:
    task.ext.when == null || task.ext.when

    script:
    def basename   = file.getName()
    def remoteName = meta.remote_name ?: basename
    def remoteDir  = remote_dir ? (remote_dir.endsWith('/') ? remote_dir : remote_dir + '/') : ''
    def remote_path = "${remoteDir}${remoteName}"
    """
    set -euo pipefail
    : "${WEBIN_USER:?WEBIN_USER not set}"
    : "${WEBIN_PASSWORD:?WEBIN_PASSWORD not set}"
    command -v ascp >/dev/null 2>&1 || { echo 'ascp is required' >&2; exit 127; }

    # Upload file and md5 to Webin drop-box via Aspera
    ascp -QT -l ${ascp_limit} -P 33001 -k 1 "${file}" "${WEBIN_USER}@${ascp_host}:/${remote_path}"
    ascp -QT -l ${ascp_limit} -P 33001 -k 1 "${md5}"  "${WEBIN_USER}@${ascp_host}:/${remote_path}.md5"
    """
}
