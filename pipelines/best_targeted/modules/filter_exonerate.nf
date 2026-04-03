// FILTER_EXONERATE
// Filter raw exonerate GFF2 output by coverage and percent identity.
// Produces per-batch GFF3.

process FILTER_EXONERATE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta), path(exonerate_gff)
    val   query_type    // 'cdna' or 'protein'
    path  query_fasta   // optional; pass [] for approximate coverage

    output:
    tuple val(meta), path("*.filtered.gff3"), emit: gff3
    path  "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: meta.id
    def args        = task.ext.args   ?: ''
    def min_cov     = params.exonerate_min_coverage  ?: 50
    def min_pid     = params.exonerate_min_pid       ?: 50
    def query_arg   = query_fasta ? "--query_fasta ${query_fasta}" : ''
    """
    filter_exonerate.py \\
        --exonerate_gff ${exonerate_gff} \\
        --out           ${prefix}.filtered.gff3 \\
        --query_type    ${query_type} \\
        --min_coverage  ${min_cov} \\
        --min_pid       ${min_pid} \\
        --sample_id     ${meta.id} \\
        ${query_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\\n' > ${prefix}.filtered.gff3
    printf 'chr1\\texonerate\\tgene\\t1000\\t5000\\t500\\t+\\t.\\tID=${prefix}_bt_gene_00000001;Name=NM_001;biotype=cdna_alignment;coverage=95.0;pid=97.0\\n' \\
        >> ${prefix}.filtered.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
