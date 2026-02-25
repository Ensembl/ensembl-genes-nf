process ENA_GENERATE_XML {
    label 'process_light'
    tag { meta.id }
    publishDir "${params.outdir}/ena_submission/xml", mode: 'copy', pattern: '*.xml', saveAs: { fn -> "${meta.id}/$fn" }

    input:
    tuple val(meta), val(row), path(file), path(md5), val(remote_path)
    val hold_until

    output:
    tuple val(meta), path('analysis.xml'), path('submission.xml'), emit: xml

    script:
    def hdr  = row.keySet().join('\t')
    def vals = row.values().collect { it ?: '' }.join('\t')
    def analysis_type = (row.containsKey('analysis_type') && row.analysis_type) ? row.analysis_type : 'READ_ALIGNMENT'
    """
    set -euo pipefail
    outdir=xml_out
    mkdir -p "$outdir"

    cat > row.tsv <<'EOF'
    ${hdr}
    ${vals}
    EOF

    pipelines/ena_submit/bin/generate_analysis_xml.py \
      --manifest-row row.tsv \
      --remote-path "${remote_path}" \
      --md5 ${md5} \
      --analysis-type ${analysis_type} \
      --hold-until ${hold_until ?: ''} \
      --outdir "$outdir" > paths.txt

    cp "$outdir/analysis.xml" .
    cp "$outdir/submission.xml" .
    """
}
