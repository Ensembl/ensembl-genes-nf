process MAKE_ORFRATER_TFAMS {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    input:
    tuple val(meta), path(bed12)
    output:
    tuple val(meta), path('tfams'), path(bed12), emit: tfams
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    mkdir -p tfams
    python \$ORFRATER_HOME/make_tfams.py --inbed ${bed12} --tfamstem tfams/tfams --force
    printf '"%s":\n    ORF-RATER: source-pinned\n    stage: make_tfams\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p tfams
    touch tfams/tfams.h5
    printf '"stub":\n    ORF-RATER: stub\n    stage: make_tfams\n' > versions.yml
    """
}
