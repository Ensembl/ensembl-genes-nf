process STRINGTIE2_COLLAPSE {
    tag "${meta.id}:${shard}:stringtie2"
    label 'process_high_memory'
    conda 'bioconda::stringtie=2.2.3'
    // Use the Bioconda/Depot SIF directly: the previously pinned Wave image
    // resolves to linux/arm64 and cannot be converted on the amd64 HPC.
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0'

    input:
    tuple val(meta), val(shard), val(resource_class), val(mapped_reads), path(bam), path(bai)

    output:
    tuple val(meta), val(shard), path('*.gtf'), emit: gtf
    tuple val(meta), val(shard), path('*.bed'), emit: bed
    path 'versions.yml', emit: versions

    script:
    def prefix = "${meta.id}.${shard}.stringtie2"
    """
    stringtie -L -p ${task.cpus} -o ${prefix}.gtf ${bam}
    gtf_to_bed12.py ${prefix}.gtf ${prefix}.bed ${meta.id}
    printf '"%s":\\n    stringtie: 2.2.3\\n' '${task.process}' > versions.yml
    """

    stub:
    def prefix = "${meta.id}.${shard}.stringtie2"
    """
    printf 'chrStub\\tstringtie\\texon\\t1\\t4\\t.\\t+\\t.\\tgene_id "${meta.id}.g1"; transcript_id "${meta.id}.t1";\\n' > ${prefix}.gtf
    printf 'chrStub\\t0\\t4\\t${meta.id}.t1\\t0\\t+\\t0\\t4\\t0\\t1\\t4,\\t0,\\n' > ${prefix}.bed
    printf '"%s":\\n    stringtie: 2.2.3-stub\\n' '${task.process}' > versions.yml
    """
}
