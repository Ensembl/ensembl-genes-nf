process RUN_RIBOTRICER {
    tag "${meta.id}:${meta.codon ?: 'all'}:${meta.shard_id ?: 'all'}"
    label 'process_single_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0'
    input:
    tuple val(meta), path(bam), path(bai), path(read_lengths), path(psite_offsets), path(orf_index), val(start_codon)
    output:
    tuple val(meta), path('raw'), emit: raw
    path 'versions.yml', emit: versions, topic: versions
    script:
    def args = task.ext.args ?: params.args_ribotricer ?: ''
    """
    mkdir -p raw
    offset_args=''
    if [ -s ${read_lengths} ] && [ -s ${psite_offsets} ]; then
        offset_args="--read_lengths \$(cat ${read_lengths}) --psite_offsets \$(cat ${psite_offsets})"
    fi
    ribotricer detect-orfs --bam ${bam} --ribotricer_index ${orf_index} \\
        --prefix raw/ribotricer \$offset_args ${args}
    test -n "\$(find raw -type f | head -1)" || { echo 'Ribotricer produced no output' >&2; exit 1; }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ribotricer: 1.5.0
        start_codon: ${start_codon}
        offsets: external
    END_VERSIONS
    """
    stub:
    """
    mkdir -p raw
    printf 'ORF_ID\tstatus\tphase_score\ttranscript_id\tchrom\tstrand\tstart_codon\nTX1_101_190_90\ttranslating\t0.9\tTX1\tchr1\t+\t${start_codon}\n' > raw/ribotricer_translating_ORFs.tsv
    printf '"stub":\n    Ribotricer: stub\n' > versions.yml
    """
}
