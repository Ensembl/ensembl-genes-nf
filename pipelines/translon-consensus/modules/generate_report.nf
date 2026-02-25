process GENERATE_HTML_REPORT {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pip_jinja2_pandas:fee727bf7c211ccf"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(consensus_results)
    val(ucsc_session_url)

    output:
    path("translon_consensus_report.html"), emit: html_report

    script:
    """
    generate_html_report.py \\
        -i . \\
        -o translon_consensus_report.html \\
        -n "Translon Consensus Analysis" \\
        -u "${ucsc_session_url}"
    """

    stub:
    """
    echo "<html><body><h1>Translon Consensus Report (Stub)</h1></body></html>" > translon_consensus_report.html
    """
}