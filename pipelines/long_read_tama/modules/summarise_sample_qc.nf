process SUMMARISE_SAMPLE_QC {
    tag "${meta.id}:sample-qc"
    label 'process_light'
    container 'community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'
    input:
    tuple val(meta), path(pre_qc), path(alignment_counts), path(alignment_identity), path(coverage), path(intron_chains)
    path summary_script
    output:
    tuple val(meta), path('sample_qc_summary.tsv'), emit: summary
    path 'versions.yml', emit: versions
    script:
    """
    python3 "${summary_script}" "${meta.id}" sample_qc_summary.tsv "${pre_qc}" "${alignment_counts}" "${alignment_identity}" "${coverage}" "${intron_chains}"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
    stub:
    """
    printf 'sample_id\\tread_count\\ttotal_alignments\\tidentity_mean\\tmean_depth\\tunique_intron_chain_count\\n${meta.id}\\t100\\t90\\t0.975000\\t12.500000\\t4\\n' > sample_qc_summary.tsv
    printf '"%s":\\n    python: 3.11.0\\n' '${task.process}' > versions.yml
    """
}
