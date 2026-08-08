/*
 * SAMPLE_COMPLETENESS
 * Compare the expected sample set with successful outputs from key stages.
 * This deliberately records missing samples rather than failing the complete
 * run because one sample-level task was ignored after retry.
 */

process SAMPLE_COMPLETENESS {
    tag "sample-completeness"
    label 'process_light'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'sample_status.tsv'

    input:
    val expected_runs
    val acquired_runs
    val aligned_runs
    val qc_runs

    output:
    path 'sample_status.tsv', emit: status

    script:
    def expected_text = expected_runs.unique().sort().collect { it.toString() }.join('\n')
    def acquired_text = acquired_runs.unique().sort().collect { it.toString() }.join('\n')
    def aligned_text = aligned_runs.unique().sort().collect { it.toString() }.join('\n')
    def qc_text = qc_runs.unique().sort().collect { it.toString() }.join('\n')
    """
    cat > expected_runs.txt <<'EOF_EXPECTED'
${expected_text}
EOF_EXPECTED
    cat > acquired_runs.txt <<'EOF_ACQUIRED'
${acquired_text}
EOF_ACQUIRED
    cat > aligned_runs.txt <<'EOF_ALIGNED'
${aligned_text}
EOF_ALIGNED
    cat > qc_runs.txt <<'EOF_QC'
${qc_text}
EOF_QC

    printf 'Run\\tAcquisition\\tAlignment\\tQC\\tOverall\\n' > sample_status.tsv
    while IFS= read -r run; do
        [ -z "\$run" ] && continue
        acquisition=FAIL
        alignment=FAIL
        qc=FAIL
        grep -Fxq "\$run" acquired_runs.txt && acquisition=PASS
        grep -Fxq "\$run" aligned_runs.txt && alignment=PASS
        grep -Fxq "\$run" qc_runs.txt && qc=PASS
        overall=FAIL
        [ "\$qc" = PASS ] && overall=PASS
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$run" "\$acquisition" "\$alignment" "\$qc" "\$overall" >> sample_status.tsv
    done < expected_runs.txt
    """

    stub:
    """
    printf 'Run\\tAcquisition\\tAlignment\\tQC\\tOverall\\n' > sample_status.tsv
    """
}
