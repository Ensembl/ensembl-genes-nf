process IRIBO_GET_CANDIDATES {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_ultra_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(gtf), path(fasta)
    output:
    tuple val(meta), path(bam), path(bai), path('raw/candidates'), path(gtf), path(fasta), emit: candidates
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_iribo ?: ''
    """
    mkdir -p raw
    iRibo --RunMode=GetCandidateORFs --Genome=${fasta} --Annotations=${gtf} --Output=raw/candidates --Threads=${task.cpus ?: 2} ${args}
    test -n "\$(find raw/candidates -type f | head -1)" || { echo 'iRibo produced no candidate output' >&2; exit 1; }
    printf '"%s":\n    iRibo: source-pinned\n    stage: candidates\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw/candidates
    printf 'candidate\n' > raw/candidates/candidate_orfs
    printf '"stub":\n    iRibo: stub\n' > versions.yml
    """
}
