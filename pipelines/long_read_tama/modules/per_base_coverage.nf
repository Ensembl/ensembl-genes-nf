process PER_BASE_COVERAGE {
    publishDir path: { "${params.outdir}/per_run/${meta.id}/qc" }, mode: 'copy', saveAs: { filename -> filename == 'versions.yml' ? null : filename }
    tag "${meta.id}:per-base-coverage"
    label 'process_medium'
    container 'https://depot.galaxyproject.org/singularity/mosdepth:0.3.6--hd299d5a_0'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('coverage_summary.tsv'), emit: summary
    path '*.per-base.bed.gz', emit: per_base
    path '*.thresholds.bed.gz', emit: thresholds
    path '*.mosdepth.global.dist.txt', emit: distribution
    path '*.mosdepth.summary.txt', emit: mosdepth_summary
    path 'versions.yml', emit: versions
    script:
    def prefix = meta.id
    def threshold_headers = params.coverage_thresholds.collect { threshold -> "breadth_${threshold}x" }.join('\\t')
    def threshold_count = params.coverage_thresholds.size()
    """
    mosdepth --by 1 --fast-mode --thresholds ${params.coverage_thresholds.join(',')} "${prefix}" "${bam}"
    mean_depth=\$(awk '\$1 != "total" && \$1 != "chrom" { bases += \$2; weighted += \$6 * \$2 } END { print (bases ? weighted / bases : 0) }' "${prefix}.mosdepth.summary.txt")
    zcat "${prefix}.thresholds.bed.gz" | awk -v sample='${meta.id}' -v mean_depth="\$mean_depth" -v threshold_count=${threshold_count} 'BEGIN { OFS="\\t"; print "sample_id","reference_bases","mean_depth","${threshold_headers}" } { length = \$3 - \$2; bases += length; for (i = 1; i <= threshold_count; i++) covered[i] += (\$(3 + i) > 0 ? length : 0) } END { printf "%s\\t%d\\t%s", sample, bases, mean_depth; for (i = 1; i <= threshold_count; i++) printf "\\t%.6f", (bases ? covered[i] / bases : 0); printf "\\n" }' > coverage_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | head -n1)
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}.per-base.bed.gz ${meta.id}.thresholds.bed.gz ${meta.id}.mosdepth.global.dist.txt ${meta.id}.mosdepth.summary.txt
    printf 'sample_id\\treference_bases\\tmean_depth\\tbreadth_1x\\tbreadth_5x\\tbreadth_10x\\n${meta.id}\\t100\\t12.500000\\t0.900000\\t0.800000\\t0.700000\\n' > coverage_summary.tsv
    printf '"%s":\\n    mosdepth: 0.3.6\\n' '${task.process}' > versions.yml
    """
}
