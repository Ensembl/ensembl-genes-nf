// PARSE_ASSEMBLY_REPORT
// Parse NCBI assembly_report.txt → seq-region synonyms TSV + metadata JSON.
// The synonyms TSV is consumed by refseq_import for coordinate remapping.

process PARSE_ASSEMBLY_REPORT {
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path assembly_report

    output:
    path "*.synonyms.tsv",       emit: synonyms
    path "*.assembly_meta.json", emit: meta
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = assembly_report.name.replaceAll(/_assembly_report\.txt$/, '')
    """
    parse_assembly_report.py \\
        --report       ${assembly_report} \\
        --synonyms_out ${prefix}.synonyms.tsv \\
        --meta_out     ${prefix}.assembly_meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = assembly_report.name.replaceAll(/_assembly_report\.txt$/, '')
    """
    printf '# RefSeq_accession\\tseq_region_name\\trole\\tlength\\n' > ${prefix}.synonyms.tsv
    printf 'NC_000001.11\\tchr1\\tassembled-molecule\\t248956422\\n' >> ${prefix}.synonyms.tsv
    printf '{"assembly_metadata": {}, "sequence_summary": {"total_sequences": 1}}' > ${prefix}.assembly_meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
