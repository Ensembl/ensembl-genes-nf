process ENA_GENERATE_XML {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/python:3.11-slim'
    publishDir "${params.outdir}/ena_submission/xml", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.id}/$fn" }

    input:
    tuple val(meta), val(row), val(file_metas), path(files), path(md5s)
    val remote_dir
    val hold_until

    output:
    tuple val(meta), path('webin_submission.xml'), emit: xml
    path('analysis.xml'), emit: analysis_xml
    path('submission.xml'), emit: submission_xml

    script:
    def escape = { v -> (v ?: '').toString().replace('\t', ' ').replace('\n', ' ') }
    def hdr = row.keySet().join('\t')
    def vals = row.values().collect { escape(it) }.join('\t')
    def analysis_type = (row.containsKey('analysis_type') && row.analysis_type) ? row.analysis_type : 'REFERENCE_ALIGNMENT'
    def remoteDir = remote_dir ? remote_dir.replaceAll('/+$', '') : ''
    def hold_until_arg = hold_until ? "--hold-until ${hold_until}" : ''
    def omit_run_refs = (params.mode == 'test' && (row.containsKey('omit_run_refs_in_test') ? row.omit_run_refs_in_test?.toString()?.toLowerCase() in ['1','true','yes'] : true))
    def omit_arg = omit_run_refs ? '--omit-run-refs' : ''
    def generator = "${moduleDir}/../bin/generate_analysis_xml.py"
    def fileRows = file_metas.withIndex().findAll { fm, idx -> !fm.is_index }.collect { fm, idx ->
        def remoteName = fm.remote_name ?: files[idx].getName()
        def remotePath = remoteDir ? "${remoteDir}/${remoteName}" : remoteName
        [
            escape(remotePath),
            escape(fm.file_type ?: row.file_type ?: ''),
            escape(md5s[idx].toString()),
            escape(fm.run_accession ?: ''),
            escape(fm.sample_accession ?: ''),
            escape(fm.experiment_accession ?: '')
        ].join('\t')
    }.join('\n')
    """
set -euo pipefail
mkdir -p xml_out

cat > row.tsv <<'EOF'
${hdr}
${vals}
EOF

cat > files.tsv <<'EOF'
remote_path	file_type	md5_path	run_accession	sample_accession	experiment_accession
${fileRows}
EOF

python3 ${generator} \
  --manifest-row row.tsv \
  --files-manifest files.tsv \
  --analysis-type ${analysis_type} \
  ${omit_arg} \
  ${hold_until_arg} \
  --outdir xml_out

cp xml_out/analysis.xml .
cp xml_out/submission.xml .
cp xml_out/webin_submission.xml .
"""
}
