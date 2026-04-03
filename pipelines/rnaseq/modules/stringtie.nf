// STRINGTIE
// Transcript assembly from a sorted BAM using StringTie2.
// Produces per-sample GTF with coverage/FPKM/TPM estimates.

process STRINGTIE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::stringtie=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0' :
        'biocontainers/stringtie:2.2.3--h43eeafb_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stringtie.gtf"), emit: gtf
    path  "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def args     = task.ext.args   ?: ''
    def strand   = meta.strandedness == 'forward'  ? '--fr' :
                   meta.strandedness == 'reverse'  ? '--rf' : ''
    """
    stringtie \\
        ${bam} \\
        -o ${prefix}.stringtie.gtf \\
        -p ${task.cpus} \\
        -m ${params.stringtie_min_length} \\
        -c ${params.stringtie_min_coverage} \\
        ${strand} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '# StringTie version 2.2.3\\n' > ${prefix}.stringtie.gtf
    printf 'chr1\\tStringTie\\ttranscript\\t1000\\t5000\\t.\\t+\\t.\\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "5.000000"; FPKM "1.234"; TPM "2.345";\\n' \\
        >> ${prefix}.stringtie.gtf
    printf 'chr1\\tStringTie\\texon\\t1000\\t2000\\t.\\t+\\t.\\tgene_id "STRG.1"; transcript_id "STRG.1.1";\\n' \\
        >> ${prefix}.stringtie.gtf
    printf 'chr1\\tStringTie\\texon\\t3000\\t5000\\t.\\t+\\t.\\tgene_id "STRG.1"; transcript_id "STRG.1.1";\\n' \\
        >> ${prefix}.stringtie.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: 2.2.3
    END_VERSIONS
    """
}
