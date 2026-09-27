process IRIBO_GENERATE_PROFILE {
    tag "${meta.id}:${meta.shard_id ?: 'all'}"
    label 'process_ultra_high'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    input:
    tuple val(meta), path(bam), path(bai), path(candidates), path(gtf), path(fasta)
    output:
    tuple val(meta), path('raw/profile'), path('raw/candidates'), emit: profile
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_iribo ?: ''
    """
    mkdir -p raw
    cp -r ${candidates} raw/candidates
    printf '%s\n' '${bam}' > raw/riboseq_bams.txt
    iRibo --RunMode=GenerateTranslationProfile --Genome=${fasta} --Annotations=${gtf} --Riboseq=raw/riboseq_bams.txt --CandidateORFs=raw/candidates/candidate_orfs --Output=raw/profile --Threads=${task.cpus ?: 2} ${args}
    test -n "\$(find raw/profile -type f | head -1)" || { echo 'iRibo produced no profile output' >&2; exit 1; }
    printf '"%s":\n    iRibo: source-pinned\n    stage: profile\n' '${task.process}' > versions.yml
    """
    stub:
    """
    mkdir -p raw/profile raw/candidates
    touch raw/profile/translation_calls raw/profile/null_distribution raw/candidates/candidate_orfs
    printf '"stub":\n    iRibo: stub\n' > versions.yml
    """
}
