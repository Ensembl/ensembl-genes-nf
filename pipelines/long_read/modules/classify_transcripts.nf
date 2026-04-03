// Classify collapsed transcript models by protein support level.
// Reads blastp hits against UniProt and updates GFF3 biotype attributes.
// Mirrors HiveClassifyTranscriptSupport with classification_type='long_read'.
//
// Biotypes assigned:
//   isoseq_supported  — blastp hit with evalue <= threshold
//   isoseq            — no blastp hit (unannotated long-read model)

process CLASSIFY_TRANSCRIPTS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    input:
    tuple val(meta),  path(gff3)
    tuple val(meta2), path(blast_tsv)

    output:
    tuple val(meta), path("*.classified.gff3"), emit: gff3
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    classify_transcripts.py \\
        --gff3      ${gff3} \\
        --blast     ${blast_tsv} \\
        --out       ${prefix}.classified.gff3 \\
        --sample-id ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.classified.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
