process PREPARE_RIBOTRICER_OFFSETS {
    tag "${meta.id}"
    label 'process_low'
    container 'python:3.12-bookworm'
    input:
    tuple val(meta), path(offsets)
    output:
    tuple val(meta), path('ribotricer_read_lengths.txt'), path('ribotricer_psite_offsets.txt'), emit: prepared
    script:
    """
    python3 ${params.translon_analysis_bin}/prepare_ribotricer_offsets.py --input ${offsets} --read-lengths ribotricer_read_lengths.txt --offsets ribotricer_psite_offsets.txt
    printf '"%s":\n    offset_adapter: python3\n' '${task.process}' > versions.yml
    """
    stub:
    """
    printf '28,30\n' > ribotricer_read_lengths.txt
    printf '13,15\n' > ribotricer_psite_offsets.txt
    printf '"stub":\n    offset_adapter: stub\n' > versions.yml
    """
}
