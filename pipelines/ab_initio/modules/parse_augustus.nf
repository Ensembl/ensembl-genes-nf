// PARSE_AUGUSTUS
// Convert raw Augustus GFF output to Ensembl-style GFF3 with ab_initio biotype.
// Filters by minimum gene length.

process PARSE_AUGUSTUS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta), path(augustus_gff)

    output:
    tuple val(meta), path("*.ab_initio.gff3"), emit: gff3
    path  "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    parse_augustus_gff.py \\
        --input          ${augustus_gff} \\
        --out            ${prefix}.ab_initio.gff3 \\
        --min_gene_length ${params.min_gene_length} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\\n' > ${prefix}.ab_initio.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
