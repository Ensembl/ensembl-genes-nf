process BED12_TO_GTF {
    tag { "${meta.pair_id}_${source_label}" }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/inputs" }, mode: 'copy', pattern: '*.gtf'

    input:
    tuple val(meta), val(source_label), path(bed12)

    output:
    tuple val(meta), val(source_label), path("${source_label}.gtf"), emit: gtf
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py bed12gtf --input '${bed12}' --output '${source_label}.gtf'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf '## synthetic BED12 conversion\\n' > '${source_label}.gtf'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
