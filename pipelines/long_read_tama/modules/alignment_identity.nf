process ALIGNMENT_IDENTITY {
    tag "${meta.id}:alignment-identity"
    label 'process_medium'
    container 'community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'
    input:
    tuple val(meta), path(bam), path(bai)
    path identity_script
    output:
    tuple val(meta), path('alignment_identity.tsv'), emit: stats
    path 'per_read_identity.tsv.gz', emit: per_read
    path 'versions.yml', emit: versions
    script:
    """
    python3 "${identity_script}" "${bam}" "${meta.id}" alignment_identity.tsv per_read_identity.tsv
    gzip -f per_read_identity.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
    stub:
    """
    printf 'sample_id\\tprimary_mapped_reads\\tidentity_weighted_mean\\tidentity_mean\\tidentity_median\\tidentity_p10\\tidentity_p90\\tidentity_min\\tidentity_max\\tidentity_ge_90\\tidentity_ge_95\\tidentity_ge_98\\tidentity_ge_99\\n${meta.id}\\t9\\t0.975000\\t0.975000\\t0.975000\\t0.950000\\t0.990000\\t0.900000\\t1.000000\\t9\\t8\\t5\\t1\\n' > alignment_identity.tsv
    printf 'read_id\\tcontig\\tidentity\\tmatches\\tsubstitutions\\tinserted_bases\\tdeleted_bases\\taligned_query_bases\\tmapq\\n' | gzip -c > per_read_identity.tsv.gz
    printf '"%s":\\n    python: 3.11.0\\n    pysam: 0.22.1\\n' '${task.process}' > versions.yml
    """
}
