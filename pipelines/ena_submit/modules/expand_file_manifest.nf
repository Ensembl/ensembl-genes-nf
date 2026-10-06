process ENA_EXPAND_FILE_MANIFEST {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/python:3.11.13-slim-bookworm'

    input:
    tuple val(meta), val(row), path(files_manifest)
    path expander_script

    output:
    tuple val(meta), val(row), path('expanded_files.tsv'), emit: expanded
    path 'versions.yml', emit: versions

    script:
    def defaultFileType = row.file_type ?: ''
    def hdr = row.keySet().join('\t')
    def vals = row.values().collect { value -> (value ?: '').toString().replace('\t', ' ').replace('\n', ' ') }.join('\t')
    """
set -euo pipefail
cat > analysis.tsv <<'EOF'
${hdr}
${vals}
EOF

python3 ${expander_script} \
  --files-tsv "${files_manifest}" \
  --analysis-tsv analysis.tsv \
  --analysis-id "${meta.id}" \
  --project-alias "${meta.project_alias}" \
  --assembly "${meta.assembly}" \
  --release "${meta.release}" \
  --default-file-type "${defaultFileType}" \
  --out expanded_files.tsv
python3 --version 2>&1 | awk '{print "ENA_EXPAND_FILE_MANIFEST:\\n  python: \"" \$2 "\""}' > versions.yml
"""

    stub:
    """
    printf '%s\n' 'analysis_id	project_alias	assembly	release	study	umbrella_study	analysis_alias	title	description	assembly_accession	reference_fasta	assembly_report	reference_supplement	last_geneset_update	partial_release_label	species	taxon_id	ref_seqs	analysis_links	analysis_attributes	analysis_type	omit_run_refs_in_test	file_path	file_type	remote_name	run_accession	sample_accession	experiment_accession' > expanded_files.tsv
    printf 'ENA_EXPAND_FILE_MANIFEST:\n  python: "stub"\n' > versions.yml
    """
}
