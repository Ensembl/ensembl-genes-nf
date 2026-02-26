process ENA_SUBMIT_WEBIN {
    label 'process_light'
    tag { meta.id }
    container 'community.wave.seqera.io/library/curl:8.18.0--78f80c4b644630b0'

    input:
    tuple val(meta), path(webin_xml)
    val webin_base
    val webin_user
    val webin_password

    output:
    tuple val(meta), path('*.queue.json'), emit: queued

    shell:
    """
    set -euo pipefail
    curl -sS -u "!{webin_user}:!{webin_password}" \
      -H 'Content-Type: application/xml' \
      --data-binary @!{webin_xml} \
      "!{webin_base}/submit/queue" > "!{meta.id}.queue.json"
    """
}
