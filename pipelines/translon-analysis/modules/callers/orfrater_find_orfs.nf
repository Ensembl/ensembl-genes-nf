process FIND_ORFRATER_ORFS {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_single_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-orfrater:1.0.0'
    input:
    tuple val(meta), path(tfams), path(bed12)
    path fasta
    output:
    tuple val(meta), path('orf.h5'), path(bed12), emit: orfs
    path 'versions.yml', emit: versions, topic: versions
    script:
    def startCodons = (params.start_codons ?: 'ATG').split(',').collect { codon -> codon.trim().toUpperCase() }.findAll { codon -> codon ==~ /[ACGT]{3}/ }
    def codonArgs = startCodons.collect { codon -> "'${codon}'" }.join(' ')
    """
    python \$ORFRATER_HOME/find_orfs_and_types.py ${fasta} --tfamstem ${tfams}/tfams --inbed ${bed12} --orfstore orf.h5 --codons ${codonArgs} --force
    test -s orf.h5
    printf '"%s":\n    ORF-RATER: source-pinned\n    stage: find_orfs\n' '${task.process}' > versions.yml
    """
    stub:
    """
    touch orf.h5
    printf '"stub":\n    ORF-RATER: stub\n    stage: find_orfs\n' > versions.yml
    """
}
