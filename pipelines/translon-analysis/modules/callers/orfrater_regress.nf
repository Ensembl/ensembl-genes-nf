process REGRESS_ORFRATER {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    input:
    tuple val(meta), path(bam), path(bed12), path(orfstore), path(offsets)
    output:
    tuple val(meta), path('regression.h5'), path('metagene.txt'), path('orf.h5'), path(bed12), emit: regression
    path 'versions.yml', emit: versions, topic: versions
    script:
    """
    cp ${offsets} offsets.txt
    python \$ORFRATER_HOME/regress_orfs.py ${bam} --orfstore ${orfstore} --inbed ${bed12} --offsetfile offsets.txt --regressfile regression.h5 --metagenefile metagene.txt --numproc ${task.cpus ?: 1} --force
    test -s regression.h5
    printf '"%s":\n    ORF-RATER: source-pinned\n    stage: regress\n' '${task.process}' > versions.yml
    """
    stub:
    """
    touch regression.h5 metagene.txt
    printf '"stub":\n    ORF-RATER: stub\n    stage: regress\n' > versions.yml
    """
}
