process ENA_GENERATE_UMBRELLA_XML {
    label 'process_light'
    tag { meta.id }
    publishDir "${params.outdir}/ena_submission/umbrella", mode: 'copy', pattern: 'umbrella_*.xml'

    input:
    tuple val(meta), path(project_accessions, stageAs: 'project_accessions.tsv'), path(analysis_accessions, stageAs: 'analysis_accessions.tsv')
    val umbrella_action
    val umbrella_alias
    val umbrella_title
    val umbrella_description

    output:
    tuple val(meta), path('umbrella_project.xml'), path('umbrella_submission.xml'), path('project_accessions.tsv'), emit: xml
    path 'versions.yml', emit: versions

    script:
    def generator = "${moduleDir}/../bin/generate_umbrella_xml.py"
    def accession_args = "--project-accessions ${project_accessions} --analysis-accessions ${analysis_accessions}"
    """
    set -euo pipefail
    python3 ${generator} \\
      --action '${umbrella_action}' \\
      --alias '${umbrella_alias}' \\
      --title '${umbrella_title}' \\
      --description '${umbrella_description}' \\
      --umbrella-accession '${params.umbrella_study}' \\
      --format dropbox \\
      ${accession_args} \\
      --outdir .
    printf 'ENA_GENERATE_UMBRELLA_XML:\\n  python: "%s"\\n' "\$(python3 --version 2>&1 | awk '{print \$2}')" > versions.yml
    """

    stub:
    """
    touch umbrella_project.xml umbrella_submission.xml
    touch project_accessions.tsv
    printf 'ENA_GENERATE_UMBRELLA_XML:\\n  python: "stub"\\n' > versions.yml
    """
}
