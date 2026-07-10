process ENA_POLL_WEBIN {
    label 'process_light'
    tag 'poll'
    container 'docker.io/library/python:3.11-slim'
    publishDir "${params.outdir}/ena_submission", mode: 'copy', pattern: 'accessions.tsv'

    input:
    val queue_responses
    val webin_user
    val webin_password
    val poll_interval
    val poll_max_attempts

    output:
    path('accessions.tsv'), emit: accessions

    script:
    def qlist = (queue_responses instanceof List) ? queue_responses : [ queue_responses ]
    def qargs = qlist.collect { "\"${it}\"" }.join(' ')
    def poller = "${moduleDir}/../bin/poll_webin.py"
    """
    set -euo pipefail
    python3 ${poller} ${qargs} \
      --webin-user "${webin_user}" \
      --webin-password "${webin_password}" \
      --interval ${poll_interval} \
      --max-attempts ${poll_max_attempts} \
      --out accessions.tsv
    """
}
