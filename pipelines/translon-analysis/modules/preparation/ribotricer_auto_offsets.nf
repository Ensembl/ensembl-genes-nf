process PREPARE_RIBOTRICER_AUTO_OFFSETS {
    tag "${meta.id}"
    label 'process_low'
    container 'python:3.12-bookworm'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path('ribotricer_read_lengths.txt'), path('ribotricer_psite_offsets.txt'), emit: prepared
    script:
    """
    : > ribotricer_read_lengths.txt
    : > ribotricer_psite_offsets.txt
    """
}
