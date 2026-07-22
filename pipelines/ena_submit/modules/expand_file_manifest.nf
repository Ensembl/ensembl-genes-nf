process ENA_EXPAND_FILE_MANIFEST {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/python:3.11-slim'

    input:
    tuple val(meta), val(row), path(files_manifest)
    path expander_script

    output:
    tuple val(meta), val(row), path('expanded_files.tsv'), emit: expanded

    script:
    def defaultFileType = row.file_type ?: ''
    def escape = { v -> (v ?: '').toString().replace('\t', ' ').replace('\n', ' ') }
    def hdr = row.keySet().join('\t')
    def vals = row.values().collect { escape(it) }.join('\t')
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
"""
}
