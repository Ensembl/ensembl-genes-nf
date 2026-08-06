process ENA_GENERATE_PROJECT_XML {
    label 'process_light'
    tag { meta.alias }
    container 'docker.io/library/python:3.11.13-slim-bookworm'
    publishDir "${params.outdir}/ena_submission/projects", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.alias}/$fn" }

    input:
    val meta

    output:
    tuple val(meta), path('webin_project.xml'), emit: xml
    path 'versions.yml', emit: versions

    script:
    def generator = "${moduleDir}/../bin/generate_project_xml.py"
    // Use shell block to avoid Groovy string interpolation issues with $ and $(...)
    """
set -euo pipefail
hold_until_arg=""
if [ -n "${meta.hold_until}" ]; then hold_until_arg="--hold-until ${meta.hold_until}"; fi
python3 ${generator} \
  --alias ${meta.alias} \
  --name "${meta.name}" \
  --title "${meta.title}" \
  --description "${meta.description}" \
  \${hold_until_arg} \
  --outdir .
printf 'ENA_GENERATE_PROJECT_XML:\n  python: "%s"\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
"""

    stub:
    """
    touch webin_project.xml
    printf 'ENA_GENERATE_PROJECT_XML:\n  python: "stub"\n' > versions.yml
    """
}
