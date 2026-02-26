process ENA_GENERATE_XML {
    label 'process_light'
    tag { meta.id }
    container 'docker.io/library/python:3.11-slim'
    publishDir "${params.outdir}/ena_submission/xml", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.id}/$fn" }

    input:
    tuple val(meta), val(row), path(file), path(md5)
    val remote_dir
    val hold_until

    output:
    tuple val(meta), path('webin_submission.xml'), emit: xml

    script:
    def hdr         = row.keySet().join('\t')
    def vals        = row.values().collect { it ?: '' }.join('\t')
    def analysis_type = (row.containsKey('analysis_type') && row.analysis_type) ? row.analysis_type : 'REFERENCE_ALIGNMENT'
    def remoteName  = meta.remote_name ?: file.getName()
    def remoteDir   = remote_dir ? remote_dir.replaceAll('/+$', '') : ''
    def remote_path   = remoteDir ? "${remoteDir}/${remoteName}" : remoteName
    def hold_until_arg = hold_until ? "--hold-until ${hold_until}" : ''
    def omit_run_refs = (params.mode == 'test' && (row.containsKey('omit_run_refs_in_test') ? row.omit_run_refs_in_test?.toString()?.toLowerCase() in ['1','true','yes'] : true))
    def omit_arg = omit_run_refs ? '--omit-run-refs' : ''
    """
    set -euo pipefail
    mkdir -p xml_out

    cat > row.tsv <<'EOF'
    ${hdr}
    ${vals}
    EOF

    generate_analysis_xml.py \
      --manifest-row row.tsv \
      --remote-path "${remote_path}" \
      --md5 ${md5} \
      --analysis-type ${analysis_type} \
      ${omit_arg} \
      ${hold_until_arg} \
      --outdir xml_out

    cp xml_out/analysis.xml .
    cp xml_out/submission.xml .
    cp xml_out/webin_submission.xml .
    """
}
