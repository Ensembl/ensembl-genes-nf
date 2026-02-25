process ENA_SUBMIT_WEBIN {
    label 'process_light'
    tag { meta.id }
    publishDir "${params.outdir}/ena_submission/receipts/${meta.id}", mode: 'copy', pattern: 'receipt.xml'

    input:
    tuple val(meta), path(analysis_xml), path(submission_xml)
    val submit_api
    val webin_v1_base
    val webin_v2_base

    output:
    tuple val(meta), path('receipt.xml'), emit: receipt

    shell:
    """
    set -euo pipefail
    : "${WEBIN_USER:?WEBIN_USER not set}"
    : "${WEBIN_PASSWORD:?WEBIN_PASSWORD not set}"
    command -v curl >/dev/null 2>&1 || { echo 'curl is required' >&2; exit 127; }

    URL_V1="${webin_v1_base%/}/submit/"
    URL_V2="${webin_v2_base%/}/analyses/submit/queue"

    if [[ "${submit_api}" == "v1" ]]; then
      curl -sS -u "$WEBIN_USER:$WEBIN_PASSWORD" \
        -F "SUBMISSION=@${submission_xml}" \
        -F "ANALYSIS=@${analysis_xml}" \
        "$URL_V1" > receipt.xml
    else
      echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > payload.xml
      echo "<WEBIN_SUBMISSION>" >> payload.xml
      echo "<SUBMISSION>" >> payload.xml
      sed -n 's/^<\?xml.*\?>$//; p' ${submission_xml} | sed '1d' >> payload.xml || true
      echo "</SUBMISSION>" >> payload.xml
      echo "<ANALYSIS>" >> payload.xml
      sed -n 's/^<\?xml.*\?>$//; p' ${analysis_xml} | sed '1d' >> payload.xml || true
      echo "</ANALYSIS>" >> payload.xml
      echo "</WEBIN_SUBMISSION>" >> payload.xml

      curl -sS -u "$WEBIN_USER:$WEBIN_PASSWORD" \
        -H 'Content-Type: application/xml' \
        --data-binary @payload.xml \
        "$URL_V2" > receipt.xml
    fi
    """
}
