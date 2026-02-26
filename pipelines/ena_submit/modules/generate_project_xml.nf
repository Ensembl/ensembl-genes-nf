process ENA_GENERATE_PROJECT_XML {
    label 'process_light'
    tag { meta.alias }
    container 'docker.io/library/python:3.11-slim'
    publishDir "${params.outdir}/ena_submission/projects", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.alias}/$fn" }

    input:
    tuple val(meta)

    output:
    tuple val(meta), path('webin_project.xml'), emit: xml

    script:
    def hold_until_arg = meta.hold_until ? "--hold-until ${meta.hold_until}" : ''
    """
    set -euo pipefail
    generate_project_xml.py \
      --alias ${meta.alias} \
      --name "${meta.name}" \
      --title "${meta.title}" \
      --description "${meta.description}" \
      ${hold_until_arg} \
      --outdir .
    """
}

