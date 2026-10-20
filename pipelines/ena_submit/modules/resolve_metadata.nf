process ENA_RESOLVE_METADATA {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/python:3.11.13-slim-bookworm'
    publishDir "${params.outdir}/ena_submission/metadata", mode: 'copy', pattern: 'resolved_*.tsv', saveAs: { fn -> "${meta.id}/${fn}" }

    input:
    tuple val(meta), val(row), path(files_manifest)
    path resolver_script

    output:
    tuple val(meta), val(row), path('resolved_files.tsv'), path('resolved_summary.tsv'), emit: resolved
    path 'versions.yml', emit: versions

    script:
    """
set -euo pipefail
python3 ${resolver_script} \\
  --input "${files_manifest}" \\
  --output resolved_files.tsv \\
  --summary-output resolved_summary.tsv \\
  --project-alias "${meta.project_alias}" \\
  --analysis-id "${meta.id}" \\
  --assembly-accession "${meta.assembly}" \\
  --release "${meta.release}"
python3 --version 2>&1 | awk '{print "ENA_RESOLVE_METADATA:\\n  python: \"" \\$2 "\""}' > versions.yml
"""

    stub:
    def stubSpecies = row.species ?: ''
    def stubTaxon = row.taxon_id ?: ''
    """
cp ${files_manifest} resolved_files.tsv
printf 'run_accessions\\tsample_accessions\\tbiosample_accessions\\texperiment_accessions\\tsource_study_accessions\\tscientific_name\\ttaxon_id\\taggregation_scope\\tmetadata_source\\tspecies\\tproject_alias\\tanalysis_id\\tassembly_accession\\trelease\\n' > resolved_summary.tsv
printf '\\t\\t\\t\\t\\t${stubSpecies}\\t${stubTaxon}\\t\\tstub\\t${stubSpecies}\\t${meta.project_alias}\\t${meta.id}\\t${meta.assembly}\\t${meta.release}\\n' >> resolved_summary.tsv
printf 'ENA_RESOLVE_METADATA:\\n  python: "stub"\\n' > versions.yml
"""
}
