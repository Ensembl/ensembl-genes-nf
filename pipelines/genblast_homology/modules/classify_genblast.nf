// CLASSIFY_GENBLAST
// Parse GenBlast GFF, apply PID/coverage filters, assign tiered biotypes
// (genblast_1..7) using the standard Ensembl classification, and emit GFF3
// with gene/transcript/exon hierarchy.

process CLASSIFY_GENBLAST {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.classified.gff3"), emit: gff3
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def min_pid  = params.min_pid      ?: 30
    def min_cov  = params.min_coverage ?: 50
    """
    classify_genblast.py \\
        --gff          ${gff} \\
        --out          ${prefix}.classified.gff3 \\
        --min-pid      ${min_pid} \\
        --min-coverage ${min_cov} \\
        --sample-id    ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\nchr1\tgenBlastG\tgene\t1000\t2000\t.\t+\t.\tID=${prefix}_gbh_gene_00000001;Name=STUB;biotype=genblast_1\n' \\
        > ${prefix}.classified.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
