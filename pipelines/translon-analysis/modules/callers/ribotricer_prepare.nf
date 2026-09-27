process PREPARE_RIBOTRICER_ORFS {
    tag "${meta.id}:${start_codon}:${meta.shard_id ?: 'all'}"
    label 'process_medium'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    container 'quay.io/biocontainers/ribotricer:1.5.0--pyhdfd78af_0'
    input:
    tuple val(meta), path(bam), path(bai), path(read_lengths), path(psite_offsets), val(start_codon)
    path gtf
    path fasta
    output:
    tuple val(meta), path(bam), path(bai), path(read_lengths), path(psite_offsets), path('raw/orfs_candidate_orfs.tsv'), val(start_codon), emit: index
    script:
    """
    mkdir -p raw
    ribotricer prepare-orfs --gtf ${gtf} --fasta ${fasta} --prefix raw/orfs --start_codons ${start_codon}
    test -s raw/orfs_candidate_orfs.tsv
    """
    stub:
    """
    mkdir -p raw
    printf 'orf_id\ttranscript_id\tstart\tstop\nORF1\tTX1\t1\t100\n' > raw/orfs_candidate_orfs.tsv
    """
}
