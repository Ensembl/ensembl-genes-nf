process RUN_RIBOCODE {
    tag "${meta.id}:${meta.codon ?: 'all'}:${meta.shard_id ?: 'all'}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribocode:1.2.15--pyhdc42f0e_1'
    input:
    tuple val(meta), path(bam), path(bai), path(transcriptome_gtf), path(transcriptome_fasta), val(start_codon)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_ribocode ?: ''
    """
    mkdir -p raw
    RiboCode_onestep -g ${transcriptome_gtf} -f ${transcriptome_fasta} -r ${bam} \\
        -l no -s ${start_codon} -o raw/ribocode -t ${task.cpus ?: 1} ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'RiboCode produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboCode: \$(RiboCode_onestep --version 2>&1 | head -n1 || true)
        start_codon: ${start_codon}
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'chrom\tstart\tend\ttranscript_id\tframe\tscore\tpval\nchr1\t100\t200\tTX1\t0\t10\t0.01\n' > raw/ribocode.tsv
    printf '"stub":\n    RiboCode: stub\n' > versions.yml
    """
}
