process ENA_SUBMIT_WEBIN {
    label 'process_light'
    tag { meta.id ?: meta.alias ?: 'submission' }
    container 'community.wave.seqera.io/library/curl:8.18.0--78f80c4b644630b0'
    secret 'ENA_WEBIN_PASSWORD'

    input:
    tuple val(meta), path(webin_xml)
    val webin_base
    val webin_user

    output:
    tuple val(meta), path('*.queue.json'), emit: queued

    shell:
    '''
    set -euo pipefail
    ID="!{meta.id ?: meta.alias ?: 'submission'}"
    curl -sS --fail-with-body -u "!{webin_user}:\$ENA_WEBIN_PASSWORD" \
      -H 'Content-Type: application/xml' \
      --data-binary @!{webin_xml} \
      -o "$ID.queue.json" \
      "!{webin_base}/submit/queue"
    test -s "$ID.queue.json" || { echo "ENA returned an empty queue response for $ID" >&2; exit 1; }
    '''
}
