process MSA_REFERENCE_SETUP {
    label 'process_medium'
    container params.msa_extract_container
    tag 'maf-reference'
    publishDir "${params.outdir}/00_references/maf", mode: 'copy', pattern: '*'

    input:
    path manifest
    tuple val(meta), path(instances)
    tuple val(genome_meta), path(genome)

    output:
    path 'maf', emit: maf_dir
    path 'maf_reference_report.jsonl', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    setup_maf_reference.py --manifest ${manifest} --instances ${instances} --genome ${genome} --outdir maf --report maf_reference_report.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        mafIndex: \$(mafIndex 2>&1 | head -n1 || echo unknown)
        verification: source-checksum-plus-materialised-sha256
    END_VERSIONS
    """

    stub:
    """
    mkdir maf
    printf '{"state":"stub"}\\n' > maf_reference_report.jsonl
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
        mafIndex: stub
        verification: source-checksum-plus-materialised-sha256
    END_VERSIONS
    """
}
