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
    touch expanded_files.tsv
    printf 'ENA_EXPAND_FILE_MANIFEST:\n  python: "stub"\n' > versions.yml
    """
}
