process RUN_RIBORF {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-riborf:1.0.0'
    input:
    tuple val(meta), path(reads_sam), path(genepred), path(bed12)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_riborf ?: ''
    """
    mkdir -p raw/predictions
    perl \$RIBORF_HOME/ribORF.pl -f ${reads_sam} -c ${genepred} -o raw/predictions ${args}
    test -n "\$(find raw/predictions -type f | head -1)" || { echo 'RibORF produced no native output' >&2; exit 1; }
    printf '"%s":\n    RibORF: 2.0\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'chrom\tstart\tend\ttranscript_id\tframe\tscore\nchr1\t100\t200\tTX1\t0\t0.9\n' > raw/riborf.tsv
    printf '"stub":\n    RibORF: stub\n' > versions.yml
    """
}
