process ENA_LINK_UMBRELLA {
    label 'process_light'
    tag { umbrella_accession }
    container 'community.wave.seqera.io/library/curl:8.18.0--78f80c4b644630b0'
    secret 'ENA_WEBIN_PASSWORD'
    publishDir "${params.outdir}/ena_submission/umbrella", mode: 'copy'

    input:
    tuple val(meta), path(umbrella_project), path(umbrella_submission), path(project_accessions)
    val umbrella_accession
    val webin_user
    val release_mode
    val release_date

    output:
    tuple val(meta), path('umbrella_status.tsv'), emit: status
    path 'umbrella_receipt.xml', emit: receipt
    path 'umbrella_current_project.xml', emit: current_project
    path 'release_receipt.xml', emit: release_receipt
    path 'versions.yml', emit: versions

    shell:
    '''
    set -euo pipefail

    child_accession="$(awk -F '\\t' 'NR > 1 && ($4 == "SUCCESS" || $4 == "EXISTS") && $3 ~ /^PRJ[A-Z]+[0-9]+$/ { print $3; exit }' !{project_accessions})"
    test -n "$child_accession" || { echo "No successful child project accession found" >&2; exit 1; }

    curl -sS --fail-with-body \
      -u "!{webin_user}:$ENA_WEBIN_PASSWORD" \
      -F "SUBMISSION=@!{umbrella_submission}" \
      -F "PROJECT=@!{umbrella_project}" \
      "!{params.mode == 'prod' ? 'https://www.ebi.ac.uk' : 'https://wwwdev.ebi.ac.uk'}/ena/submit/drop-box/submit/" \
      -o umbrella_receipt.xml
    grep -q 'success="true"' umbrella_receipt.xml || { cat umbrella_receipt.xml >&2; exit 1; }

    curl -sS --fail-with-body \
      -u "!{webin_user}:$ENA_WEBIN_PASSWORD" \
      "!{params.mode == 'prod' ? 'https://www.ebi.ac.uk' : 'https://wwwdev.ebi.ac.uk'}/ena/submit/drop-box/projects/!{umbrella_accession}" \
      -o umbrella_current_project.xml
    grep -q "CHILD_PROJECT accession=\\\"$child_accession\\\"" umbrella_current_project.xml || {
      echo "ENA accepted the umbrella update but the child link was not visible in the retrieved parent project" >&2
      cat umbrella_current_project.xml >&2
      exit 1
    }

    release_status="NOT_REQUESTED"
    if [ "!{release_mode}" = "immediate" ]; then
      cat > release_submission.xml <<EOF
    <SUBMISSION><ACTIONS><ACTION><RELEASE target="!{umbrella_accession}"/></ACTION><ACTION><RELEASE target="$child_accession"/></ACTION></ACTIONS></SUBMISSION>
    EOF
      curl -sS --fail-with-body \
        -u "!{webin_user}:$ENA_WEBIN_PASSWORD" \
        -F "SUBMISSION=@release_submission.xml" \
        "!{params.mode == 'prod' ? 'https://www.ebi.ac.uk' : 'https://wwwdev.ebi.ac.uk'}/ena/submit/drop-box/submit/" \
        -o release_receipt.xml
      grep -q 'success="true"' release_receipt.xml || { cat release_receipt.xml >&2; exit 1; }
      release_status="IMMEDIATE"
    elif [ "!{release_mode}" = "date" ]; then
      cat > release_submission.xml <<EOF
    <SUBMISSION><ACTIONS><ACTION><HOLD target="!{umbrella_accession}" HoldUntilDate="!{release_date}"/></ACTION><ACTION><HOLD target="$child_accession" HoldUntilDate="!{release_date}"/></ACTION></ACTIONS></SUBMISSION>
    EOF
      curl -sS --fail-with-body \
        -u "!{webin_user}:$ENA_WEBIN_PASSWORD" \
        -F "SUBMISSION=@release_submission.xml" \
        "!{params.mode == 'prod' ? 'https://www.ebi.ac.uk' : 'https://wwwdev.ebi.ac.uk'}/ena/submit/drop-box/submit/" \
        -o release_receipt.xml
      grep -q 'success="true"' release_receipt.xml || { cat release_receipt.xml >&2; exit 1; }
      release_status="DATE"
    else
      printf '%s\n' '<RELEASE status="NOT_REQUESTED"/>' > release_receipt.xml
    fi

    printf 'umbrella_accession\\tchild_accession\\trelease_mode\\tstatus\\n%s\\t%s\\t%s\\tSUCCESS\\n' \
      "!{umbrella_accession}" "$child_accession" "$release_status" > umbrella_status.tsv
    printf 'ENA_LINK_UMBRELLA:\\n  curl: "%s"\\n' "$(curl --version | awk 'NR==1 {print $2}')" > versions.yml
    '''

    stub:
    '''
    printf '<RECEIPT success="true"><PROJECT accession="!{umbrella_accession}"/></RECEIPT>\n' > umbrella_receipt.xml
    printf '<PROJECT_SET><PROJECT accession="!{umbrella_accession}"><RELATED_PROJECTS><RELATED_PROJECT><CHILD_PROJECT accession="PRJEB000001"/></RELATED_PROJECT></RELATED_PROJECTS></PROJECT></PROJECT_SET>\n' > umbrella_current_project.xml
    printf '<RELEASE status="NOT_REQUESTED"/>\n' > release_receipt.xml
    printf 'umbrella_accession\\tchild_accession\\trelease_mode\\tstatus\\n!{umbrella_accession}\\tPRJEB000001\\t!{release_mode}\\tSUCCESS\\n' > umbrella_status.tsv
    printf 'ENA_LINK_UMBRELLA:\\n  curl: "stub"\\n' > versions.yml
    '''
}
