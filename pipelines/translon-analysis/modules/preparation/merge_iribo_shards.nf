process MERGE_IRIBO_SHARDS {
    tag "${meta.id}:global"
    label 'process_high'
    container 'ghcr.io/jackcurragh/translon-iribo:1.0.0'
    input:
    tuple val(meta), path(profiles, arity: '1..*'), path(candidates, arity: '1..*')
    output:
    tuple val(meta), path('merged/profile'), path('merged/candidates'), emit: merged
    script:
    """
    mkdir -p merged/profile merged/candidates
    first=1
    for profile in ${profiles}; do
        if [ -f "\$profile/translation_calls" ]; then
            if [ \$first -eq 1 ]; then cp "\$profile/translation_calls" merged/profile/translation_calls; first=0; else tail -n +2 "\$profile/translation_calls" >> merged/profile/translation_calls; fi
        fi
        if [ -f "\$profile/null_distribution" ]; then
            if [ ! -s merged/profile/null_distribution ]; then cp "\$profile/null_distribution" merged/profile/null_distribution; else tail -n +2 "\$profile/null_distribution" >> merged/profile/null_distribution; fi
        fi
    done
    first=1
    for candidate in ${candidates}; do
        if [ -f "\$candidate/candidate_orfs" ]; then
            if [ \$first -eq 1 ]; then cp "\$candidate/candidate_orfs" merged/candidates/candidate_orfs; first=0; else tail -n +2 "\$candidate/candidate_orfs" >> merged/candidates/candidate_orfs; fi
        fi
    done
    test -s merged/profile/translation_calls && test -s merged/profile/null_distribution && test -s merged/candidates/candidate_orfs
    """
    stub:
    """
    mkdir -p merged/profile merged/candidates
    printf 'index frame0 frame_sum\n1 1 1\n' > merged/profile/translation_calls
    printf 'index scrambled0 scrambled_sum0\n1 1 1\n' > merged/profile/null_distribution
    printf 'index Gene_ID contig_str gene_intersect\n1 X chr1 X\n' > merged/candidates/candidate_orfs
    """
}
