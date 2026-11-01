process GFFCOMPARE {
    tag { "${meta.pair_id}_${direction}" }
    label 'process_low'
    container 'quay.io/biocontainers/gffcompare:0.12.6--h4ac6f70_2'
    publishDir { "${params.outdir}/qc/pairwise_annotation/${meta.pair_id}/gffcompare" }, mode: 'copy', pattern: '*'

    input:
    tuple val(meta), val(direction), path(reference_gtf), path(query_gtf)

    output:
    tuple val(meta), val(direction), path("${direction}.tmap"), emit: tmap
    tuple val(meta), val(direction), path("${direction}.stats"), emit: stats
    path 'versions.yml', emit: versions

    script:
    """
    set -euo pipefail
    gffcompare -r '${reference_gtf}' -o '${direction}' '${query_gtf}'
    cp '${direction}.annotated.gtf.tmap' '${direction}.tmap'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        gffcompare: \$(gffcompare --version 2>&1 | sed 's/gffcompare v//g')
    END_VERSIONS
    """

    stub:
    """
    printf 'qry_id\\tclass_code\\n' > '${direction}.tmap'
    printf 'stub\\n' > '${direction}.stats'
    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        gffcompare: 0.12.6
    END_VERSIONS
    """
}
