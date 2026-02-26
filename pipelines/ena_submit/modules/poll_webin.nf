process ENA_POLL_WEBIN {
    label 'process_light'
    // Multi-arch curl image
    container 'curlimages/curl:8.7.1'
    publishDir "${params.outdir}/ena_submission", mode: 'copy', pattern: 'accessions.tsv'

    input:
    path(queue_responses)   // collected list of *.queue.json files
    val webin_user
    val webin_password
    val poll_interval       // seconds between attempts
    val poll_max_attempts   // max attempts per submission

    output:
    path('accessions.tsv'), emit: accessions

    shell:
    """
    set -euo pipefail
    printf 'analysis_id\tsubmission_id\taccession\tstatus\n' > accessions.tsv
    failed=0

    for qfile in !{queue_responses}; do
        analysis_id=\${qfile%.queue.json}
        # Robustly extract poll URL even if JSON spans multiple lines
        json_flat=\$(tr -d '\n' < "\$qfile")
        poll_url=\$(printf '%s' "\$json_flat" | sed -n 's/.*"href":"\([^"]*\)".*/\1/p')
        if [ -z "\$poll_url" ]; then
            echo "ERROR: could not find poll href in \$qfile" >&2
            cat "\$qfile" >&2
            exit 3
        fi
        case "\$poll_url" in
          http://*|https://*) ;;
          *) echo "ERROR: invalid poll URL extracted: '\$poll_url' from \$qfile" >&2; exit 3;;
        esac
        submission_id=\$(grep -o '"submissionId":"[^"]*"' "\$qfile" | cut -d'"' -f4)

        echo "Polling \$submission_id (\$analysis_id) ..." >&2

        receipt=""
        for i in \$(seq 1 !{poll_max_attempts}); do
            sleep !{poll_interval}
            response=\$(curl -fSs -u "!{webin_user}:!{webin_password}" "\$poll_url" || true)
            if echo "\$response" | grep -q '<RECEIPT'; then
                receipt="\$response"
                break
            fi
            echo "  attempt \$i: pending" >&2
        done

        if [ -z "\$receipt" ]; then
            echo "ERROR: timed out polling \$submission_id" >&2
            failed=1
            printf '%s\t%s\t\tTIMEOUT\n' "\$analysis_id" "\$submission_id" >> accessions.tsv
            continue
        fi

        if echo "\$receipt" | grep -q 'success="false"'; then
            echo "ERROR: \$submission_id rejected by ENA:" >&2
            echo "\$receipt" >&2
            failed=1
            printf '%s\t%s\t\tFAILED\n' "\$analysis_id" "\$submission_id" >> accessions.tsv
        else
            accession=\$(echo "\$receipt" | grep -o 'accession="[^"]*"' | head -1 | cut -d'"' -f2)
            echo "OK: \$submission_id -> \$accession" >&2
            printf '%s\t%s\t%s\tSUCCESS\n' "\$analysis_id" "\$submission_id" "\$accession" >> accessions.tsv
        fi
    done

    exit \$failed
    """
}
