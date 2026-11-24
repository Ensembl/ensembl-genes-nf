process ENA_FTP_UPLOAD {
    label 'process_single_long'
    tag { file_meta.remote_name ?: file.getName() }
    container 'docker.io/library/python:3.11.13-bookworm'
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
    python3 - "${file}" "${md5}" "${ftp_host}" "${webin_user}" \
        "${remote_path}" "${remote_path}.md5" <<'PY'
import ftplib
import os
import sys
import time

local_file, local_md5, host, user, remote_file, remote_md5 = sys.argv[1:]
timeout = 1800


def upload(ftp, home, local_path, remote_path):
    ftp.cwd(home)
    parts = remote_path.strip('/').split('/')
    remote_name = parts.pop()
    for directory in parts:
        try:
            ftp.cwd(directory)
        except ftplib.error_perm:
            try:
                ftp.mkd(directory)
            except ftplib.error_perm:
                pass
            ftp.cwd(directory)
    try:
        ftp.voidcmd('TYPE I')
        if ftp.size(remote_name) == os.path.getsize(local_path):
            print(f'Skipping existing complete file: {remote_path}', file=sys.stderr)
            return
    except ftplib.error_perm:
        pass
    with open(local_path, 'rb') as handle:
        ftp.storbinary(f'STOR {remote_name}', handle, blocksize=8 * 1024 * 1024)


def upload_all():
    with ftplib.FTP() as ftp:
        ftp.connect(host, 21, timeout=timeout)
        ftp.login(user, os.environ['ENA_WEBIN_PASSWORD'])
        ftp.set_pasv(True)
        ftp.sock.settimeout(timeout)
        home = ftp.pwd()
        upload(ftp, home, local_file, remote_file)
        upload(ftp, home, local_md5, remote_md5)


retryable = (ConnectionRefusedError, ConnectionResetError, TimeoutError,
             BrokenPipeError, ftplib.error_temp)
for attempt in range(1, 6):
    try:
        upload_all()
        break
    except retryable as exc:
        if attempt == 5:
            raise
        delay = 30 * (2 ** (attempt - 1))
        print(f'FTP attempt {attempt} failed: {exc}; retrying in {delay}s', file=sys.stderr)
        time.sleep(delay)
PY
    printf 'ENA_FTP_UPLOAD:\n  python: "%s"\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    touch ${md5}
    printf 'ENA_FTP_UPLOAD:\n  python: "stub"\n' > versions.yml
    """
}
