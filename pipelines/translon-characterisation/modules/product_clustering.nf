process PRODUCT_CLUSTERING {
    label 'process_high'
    // Production config must supply an immutable image that contains both
    // MMseqs2 and Foldseek; a score is never fabricated when structures are absent.
    container params.product_clustering_container
    tag "${meta.id}"
    publishDir "${params.outdir}/08_product_clustering", mode: 'copy', pattern: '*.clusters.jsonl'

    input:
    tuple val(meta), path(products)

    output:
    tuple val(meta), path('*.clusters.jsonl'), emit: clusters
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    product_clustering.py --instances ${products} --output ${prefix}.clusters.jsonl --mmseqs mmseqs --foldseek foldseek --min-seq-id ${params.product_cluster_min_seq_id}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version 2>/dev/null || echo unknown)
        foldseek: \$(foldseek version 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.clusters.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: 15-6f452
        foldseek: stub
    END_VERSIONS
    """
}
