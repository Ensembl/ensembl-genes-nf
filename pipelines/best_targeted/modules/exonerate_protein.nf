// EXONERATE_PROTEIN
// Align a batch of protein sequences against the softmasked genome using
// the protein2genome model.

process EXONERATE_PROTEIN {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::exonerate=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exonerate:2.4.0--hd03093a_6' :
        'biocontainers/exonerate:2.4.0--hd03093a_6' }"

    input:
    tuple val(meta), path(protein_fasta)
    path  genome_fasta

    output:
    tuple val(meta), path("*.protein.gff"), emit: gff
    path  "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: meta.id
    def args    = task.ext.args   ?: ''
    def options = '--model protein2genome --forwardcoordinates FALSE ' +
                  '--softmasktarget TRUE --exhaustive FALSE --bestn 1 ' +
                  '--showtargetgff'
    """
    exonerate \\
        ${options} \\
        ${args} \\
        --query ${protein_fasta} \\
        --target ${genome_fasta} \\
        > ${prefix}.protein.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exonerate: \$(exonerate --version 2>&1 | head -1 | sed 's/exonerate version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 2\\n' > ${prefix}.protein.gff
    printf 'chr1\\texonerate:protein2genome\\tgene\\t1200\\t4800\\t300\\t+\\t.\\tsequence P12345 ; score 300\\n' \\
        >> ${prefix}.protein.gff
    printf 'chr1\\texonerate:protein2genome\\texon\\t1200\\t2200\\t.\\t+\\t.\\tsequence P12345\\n' \\
        >> ${prefix}.protein.gff
    printf 'chr1\\texonerate:protein2genome\\texon\\t3000\\t4800\\t.\\t+\\t.\\tsequence P12345\\n' \\
        >> ${prefix}.protein.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exonerate: 2.4.0
    END_VERSIONS
    """
}
