process RUN_ORFRATER {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(bed12), path(model)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_orfrater ?: ''
    """
    mkdir -p raw
    test -n "\$(find ${model} -type f | head -1)" || { echo 'ORF-RATER model bundle is empty' >&2; exit 2; }
    python \$ORFRATER_HOME/quantify_orfs.py ${bam} --inbed ${bed12} --subdir raw --ratingsfile ${model}/orfratings.h5 --metagenefile ${model}/metagene.txt --offsetfile ${model}/offsets.txt ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'ORF-RATER produced no native output' >&2; exit 1; }
    printf '"%s":\n    ORF-RATER: source-pinned\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw
    printf 'chrom\tstart\tend\ttranscript_id\tframe\tscore\nchr1\t100\t200\tTX1\t0\t0.9\n' > raw/orfrater.tsv
    printf '"stub":\n    ORF-RATER: stub\n' > versions.yml
    """
}
