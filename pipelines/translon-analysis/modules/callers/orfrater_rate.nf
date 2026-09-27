process RATE_ORFRATER {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    input:
    tuple val(meta), path(regression), path(metagene), path(orfstore), path(bed12), path(offsets)
    output:
    tuple val(meta), path('model'), emit: model
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    mkdir -p model
    cp ${metagene} model/metagene.txt
    cp ${offsets} model/offsets.txt
    python \$ORFRATER_HOME/rate_regression_output.py ${regression} --goldallcodons --orfstore ${orfstore} --ratingsfile model/orfratings.h5 --numproc ${task.cpus ?: 1} --force
    test -s model/orfratings.h5
    printf '"%s":\n    ORF-RATER: source-pinned\n    stage: rate\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p model
    touch model/orfratings.h5 model/metagene.txt model/offsets.txt
    printf '"stub":\n    ORF-RATER: stub\n    stage: rate\n' > versions.yml
    """
}
