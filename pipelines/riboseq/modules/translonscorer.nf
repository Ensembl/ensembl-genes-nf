/*
 * TRANSLONSCORER MODULE
 * Runs TranslonScorer on per-sample BigWigs to score ORFs and produce an HTML report.
 * Inputs expect stranded BigWigs from POST_PROCESSING (forward/reverse) for a chosen bam_type.
 */

process TRANSLONSCORER {
    tag "${meta.id}"
    label 'process_medium'

    // Prefer pinned container; falls back to pip install if not present
    container "${ params.translonscorer_container ?: 'ghcr.io/jackcurragh/translonscorer:latest' }"

    input:
    // Either a single BigWig or a list [forward.bw, reverse.bw]
    tuple val(meta), path(bigwigs)
    path gtf
    path fasta

    output:
    tuple val(meta), path("${meta.id}_orfs_scored.csv"), emit: scored_csv
    tuple val(meta), path("${meta.id}_report.html"),   emit: report_html

    when:
    // Allow toggling from params
    (params.run_translonscorer ?: false)

    script:
    // User-tunable options
    def plotrng = params.translonscorer_plot_range ?: 30

    // Resolve BigWig arguments
    // bigwigs can be a single path or list; prefer stranded if two files are provided
    def fw = null
    def rv = null
    if (bigwigs instanceof List && bigwigs.size() >= 2) {
        // Try to detect by filename
        fw = bigwigs.find { value -> value.toString().toLowerCase().contains('forward') } ?: bigwigs[0]
        rv = bigwigs.find { value -> value.toString().toLowerCase().contains('reverse') } ?: bigwigs[1]
    }

    def out_base = meta.id

    // Install TranslonScorer only if not found (useful when running without container)
    def pip_spec = params.translonscorer_pip_spec ?: "TranslonScorer[full]@git+https://github.com/JackCurragh/TranslonScorer"

    if (fw && rv) {
        """
        if ! command -v translonscorer >/dev/null 2>&1; then
          python -m pip install -q --no-cache-dir ${pip_spec}
        fi
        translonscorer all \
          --forward-bigwig ${fw} \
          --reverse-bigwig ${rv} \
          --stranded \
          -s ${fasta} \
          -a ${gtf} \
          --plot-range ${plotrng} \
          --output ${out_base} \
          --log-level INFO
        """
    } else {
        def bw = bigwigs instanceof List ? bigwigs[0] : bigwigs
        """
        if ! command -v translonscorer >/dev/null 2>&1; then
          python -m pip install -q --no-cache-dir ${pip_spec}
        fi
        translonscorer all \
          -w ${bw} \
          -s ${fasta} \
          -a ${gtf} \
          --plot-range ${plotrng} \
          --output ${out_base} \
          --log-level INFO
        """
    }

    stub:
    """
    echo "sample_id,orf_id,chrom,start,end,strand,score" > ${meta.id}_orfs_scored.csv
    echo "<html><body><h1>${meta.id} TranslonScorer Report (stub)</h1></body></html>" > ${meta.id}_report.html
    """
}
