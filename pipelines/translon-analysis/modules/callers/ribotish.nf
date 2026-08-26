process RUN_RIBOTISH {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribotish:0.2.8--pyhdfd78af_0'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    path fasta
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_ribotish ?: ''
    """
    mkdir -p raw
    ribotish predict -b ${bam} -g ${gtf} -f ${fasta} -o raw/ribotish.txt --blocks --seq --aaseq -p ${task.cpus ?: 1} ${args}
    test -s raw/ribotish.txt || { echo 'Ribo-TISH produced no output' >&2; exit 1; }
    printf '"%s":\n    RiboTISH: 0.2.8\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'Tid\tGenomePos\tRiboPStatus\nTX1\tchr1:101-220:+\tT\n' > raw/ribotish.txt
    printf '"stub":\n    RiboTISH: stub\n' > versions.yml
    """
}
