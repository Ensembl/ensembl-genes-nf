process PREPARE_GFFCOMPARE_INPUTS {
    tag { meta.pair_id }
    label 'process_low'
    container 'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/gffcompare" }, mode: 'copy', pattern: '*.gtf'

    input:
    tuple val(meta), path(source_a_models), path(source_b_models)

    output:
    tuple val(meta), path('source_a.gtf'), path('source_b.gtf'), emit: gtfs
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    pairwise_qc.py models2gtf --input '${source_a_models}' --output source_a.gtf
    pairwise_qc.py models2gtf --input '${source_b_models}' --output source_b.gtf
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch source_a.gtf source_b.gtf
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        pairwise_qc: 1.0.0
    END_VERSIONS
    """
}
