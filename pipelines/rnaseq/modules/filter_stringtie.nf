// FILTER_STRINGTIE
// Filter StringTie2 GTF by coverage/length/exon count and convert to GFF3.

process FILTER_STRINGTIE {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.rnaseq.gff3"), emit: gff3
    path  "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def min_cov  = params.stringtie_min_coverage ?: 2.0
    def min_len  = params.rnaseq_min_length       ?: 200
    def min_exon = params.rnaseq_min_exons        ?: 1
    """
    filter_stringtie.py \\
        --gtf          ${gtf} \\
        --out          ${prefix}.rnaseq.gff3 \\
        --sample       ${meta.id} \\
        --biotype      rnaseq_tissue \\
        --min_coverage ${min_cov} \\
        --min_length   ${min_len} \\
        --min_exons    ${min_exon} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 3\\n' > ${prefix}.rnaseq.gff3
    printf 'chr1\\tStringTie2\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=${prefix}_rna_gene_STRG.1;biotype=rnaseq_tissue\\n' \\
        >> ${prefix}.rnaseq.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
