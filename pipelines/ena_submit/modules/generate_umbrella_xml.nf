process ENA_GENERATE_UMBRELLA_XML {
    label 'process_light'
    tag { meta.alias }
    container 'docker.io/library/python:3.11.13-slim-bookworm'
    publishDir "${params.outdir}/ena_submission/umbrella", mode: 'copy', pattern: 'umbrella.xml'

    input:
    tuple val(meta), path(project_accessions)
    val umbrella_accession
    val umbrella_title
    val umbrella_description

    output:
    tuple val(meta), path('umbrella.xml'), emit: xml
    path 'versions.yml', emit: versions

    script:
    def generator = "${moduleDir}/../bin/generate_umbrella_xml.py"
    """
set -euo pipefail
python3 ${generator} \\
  --umbrella-accession "${umbrella_accession}" \\
  --child-accessions "${project_accessions}" \\
  --title "${umbrella_title}" \\
  --description "${umbrella_description}" \\
  --outdir .
printf 'ENA_GENERATE_UMBRELLA_XML:\\n  python: "%s"\\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
"""

    stub:
    """
touch umbrella.xml
printf 'ENA_GENERATE_UMBRELLA_XML:\\n  python: "stub"\\n' > versions.yml
"""
}
