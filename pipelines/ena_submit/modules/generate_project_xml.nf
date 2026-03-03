process ENA_GENERATE_PROJECT_XML {
    label 'process_light'
    tag { meta.alias }
    container 'docker.io/library/python:3.11-slim'
    publishDir "${params.outdir}/ena_submission/projects", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.alias}/$fn" }

    input:
    val meta

    output:
    tuple val(meta), path('webin_project.xml'), emit: xml

    // Use shell block to avoid Groovy string interpolation issues with $ and $(...)
    shell:
    '''
    set -euo pipefail
    hold_until_arg=""
    if [ -n "!{meta.hold_until}" ]; then hold_until_arg="--hold-until !{meta.hold_until}"; fi
    python3 $(command -v generate_project_xml.py) \
      --alias !{meta.alias} \
      --name "!{meta.name}" \
      --title "!{meta.title}" \
      --description "!{meta.description}" \
      ${hold_until_arg} \
      --outdir .
    '''
}
