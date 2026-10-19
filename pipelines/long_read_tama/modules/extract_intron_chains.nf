process EXTRACT_INTRON_CHAINS {
    tag "${meta.id}:intron-chains"
    label 'process_medium'
    container 'community.wave.seqera.io/library/flair_coreutils_python:9744ce1f6049d5bd'
    input:
    tuple val(meta), path(bam), path(bai)
    path chain_script
    output:
    tuple val(meta), path('intron_chain_summary.tsv'), emit: summary
    path 'intron_chains.tsv.gz', emit: chains
    path 'junctions.bed', emit: junctions
    path 'versions.yml', emit: versions
    script:
    """
    python3 "${chain_script}" "${bam}" "${meta.id}" intron_chain_summary.tsv intron_chains.tsv junctions.bed
    gzip -f intron_chains.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
    stub:
    """
    printf 'sample_id\\tprimary_mapped_reads\\tspliced_reads\\tunique_intron_chain_count\\tunique_junction_count\\tmean_introns_per_spliced_read\\n${meta.id}\\t9\\t7\\t4\\t8\\t2.000000\\n' > intron_chain_summary.tsv
    printf 'read_id\\tcontig\\tstrand\\tintron_count\\tintron_chain\\tmapq\\n' | gzip -c > intron_chains.tsv.gz
    touch junctions.bed
    printf '"%s":\\n    python: 3.11.0\\n    pysam: 0.22.1\\n' '${task.process}' > versions.yml
    """
}
