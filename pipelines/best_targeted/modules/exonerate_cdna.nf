// EXONERATE_CDNA
// Align a batch of cDNA sequences against the softmasked genome using
// the cdna2genome model.  Input genome is SOFTMASKED.

process EXONERATE_CDNA {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::exonerate=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exonerate:2.4.0--hd03093a_6' :
        'biocontainers/exonerate:2.4.0--hd03093a_6' }"

    input:
    tuple val(meta), path(cdna_fasta)
    path  genome_fasta

    output:
    tuple val(meta), path("*.cdna.gff"), emit: gff
    path  "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: meta.id
    def args    = task.ext.args   ?: ''
    def options = '--model cdna2genome --forwardcoordinates FALSE ' +
                  '--softmasktarget TRUE --exhaustive FALSE ' +
                  '--score 500 --saturatethreshold 100 ' +
                  '--dnawordlen 15 --codonwordlen 15 ' +
                  '--dnahspthreshold 60 --bestn 10 ' +
                  '--showtargetgff'
    """
    exonerate \\
        ${options} \\
        ${args} \\
        --query ${cdna_fasta} \\
        --target ${genome_fasta} \\
        > ${prefix}.cdna.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exonerate: \$(exonerate --version 2>&1 | head -1 | sed 's/exonerate version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##gff-version 2\\n' > ${prefix}.cdna.gff
    printf 'chr1\\texonerate:cdna2genome\\tgene\\t1000\\t5000\\t500\\t+\\t.\\tsequence NM_001 ; score 500\\n' \\
        >> ${prefix}.cdna.gff
    printf 'chr1\\texonerate:cdna2genome\\texon\\t1000\\t2000\\t.\\t+\\t.\\tsequence NM_001\\n' \\
        >> ${prefix}.cdna.gff
    printf 'chr1\\texonerate:cdna2genome\\texon\\t3000\\t5000\\t.\\t+\\t.\\tsequence NM_001\\n' \\
        >> ${prefix}.cdna.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        exonerate: 2.4.0
    END_VERSIONS
    """
}
