process TMERGE_COLLAPSE {
    tag "${meta.id}:${shard}:tmerge"
    label 'process_high_memory'
    container 'community.wave.seqera.io/library/pip_pyfaidx_setuptools_six_pruned:b4fcc6bbecfa2bea'

    input:
    tuple val(meta), val(shard), path(reads_gtf)

    output:
    tuple val(meta), val(shard), path('*.gtf'), emit: gtf
    tuple val(meta), val(shard), path('*.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    tmerge ${params.tmerge_args} --input ${reads_gtf} --output ${prefix}.gtf
    gtf_to_bed12.py ${prefix}.gtf ${prefix}.bed ${meta.id}
    printf '"%s":\\n    tmerge: runtime\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = "${meta.id}.${shard}.tmerge"
    """
    printf 'chrStub\\ttmerge\\texon\\t1\\t4\\t.\\t+\\t.\\ttranscript_id "${meta.id}.t1";\\n' > ${prefix}.gtf
    printf 'chrStub\\t0\\t4\\t${meta.id}.t1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${prefix}.bed
    printf '"%s":\\n    tmerge: stub\\n' '${task.process}' > versions.yml
    """
}
