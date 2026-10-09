process ENA_POLL_WEBIN {
    label 'process_light'
    tag 'poll'
    container 'docker.io/library/python:3.11.13-slim-bookworm'
    secret 'ENA_WEBIN_PASSWORD'
    publishDir "${params.outdir}/ena_submission", mode: 'copy', pattern: 'accessions.tsv'

    input:
    val queue_responses
    val webin_user
    val poll_interval
    val poll_max_attempts

    output:
    path('accessions.tsv'), emit: accessions
    path 'versions.yml', emit: versions

    script:
    def qlist = (queue_responses instanceof List) ? queue_responses : [ queue_responses ]
    def qargs = qlist.collect { "\"${it}\"" }.join(' ')
    def poller = "${moduleDir}/../bin/poll_webin.py"
    """
    set -euo pipefail
    python3 ${poller} ${qargs} \
      --webin-user "${webin_user}" \
      --webin-password-env ENA_WEBIN_PASSWORD \
      --interval ${poll_interval} \
      --max-attempts ${poll_max_attempts} \
      --out accessions.tsv
    printf 'ENA_POLL_WEBIN:\n  python: "%s"\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    touch accessions.tsv
    printf 'ENA_POLL_WEBIN:\n  python: "stub"\n' > versions.yml
    """
}
