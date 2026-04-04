// PROJECT_TRANSCRIPTS
// Map source (reference) gene models onto a target genome via a UCSC chain file.
// Filters projected transcripts by coverage threshold.

process PROJECT_TRANSCRIPTS {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${params.outdir}/projection", mode: 'copy', overwrite: true,
               pattern: '*.projected.gff3'

    input:
    path chain
    path source_gff3

    output:
    path "projected.gff3",  emit: gff3
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    project_transcripts.py \\
        --source_gff3  ${source_gff3} \\
        --chain        ${chain} \\
        --out          projected.gff3 \\
        --min_coverage ${params.min_coverage} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > projected.gff3
    printf 'chr1\\tprojection\\tgene\\t1000\\t5000\\t.\\t+\\t.\\tID=proj_gene_00000001;biotype=projected_transcript\\n' \\
        >> projected.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
