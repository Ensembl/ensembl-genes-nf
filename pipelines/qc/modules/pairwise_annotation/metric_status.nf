process METRIC_STATUS {
    tag { meta.pair_id }
    label 'process_light'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/summary" }, mode: 'copy', pattern: 'metric_status.tsv'

    input:
    tuple val(meta), val(status_entries)

    output:
    tuple val(meta), path('metric_status.tsv'), emit: status
    path 'versions.yml', emit: versions

    script:
    def args = status_entries.collect { metric, status, reason -> "--entry '${metric}' '${status}' '${reason}'" }.join(' ')
    """
    set -euo pipefail
    pairwise_qc.py status --pair-id '${meta.pair_id}' --output metric_status.tsv ${args}
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    def stubRows = status_entries.collect { metric, status, reason -> "${meta.pair_id}\\t${metric}\\t${status}\\t${reason}" }.join('\\n')
    """
    printf 'pair_id\\tmetric\\tstatus\\treason\\n' > metric_status.tsv
    printf '${stubRows}\\n' >> metric_status.tsv
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
