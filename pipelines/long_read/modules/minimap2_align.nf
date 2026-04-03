process MINIMAP2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::minimap2=2.28 bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_3' :
        'biocontainers/minimap2:2.28--he4a0461_3' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(reference)
    val   bam_format
    val   bam_index_extension
    val   cigar_paf_format
    val   cigar_bam

    output:
    tuple val(meta), path("*.paf"),     optional: true, emit: paf
    tuple val(meta), path("*.bam"),     optional: true, emit: bam
    tuple val(meta), path("*.bam.*"),   optional: true, emit: index
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: meta.id
    def args       = task.ext.args   ?: ''
    def bam_output = bam_format       ? "| samtools view -bS -" : ""
    def cs_flag    = cigar_paf_format ? "--cs=long" : ""
    def L_flag     = cigar_bam        ? "-L" : ""
    def a_flag     = bam_format       ? "-a" : ""
    """
    minimap2 \\
        ${args} \\
        ${a_flag} \\
        ${cs_flag} \\
        ${L_flag} \\
        -t ${task.cpus} \\
        ${reference} \\
        ${reads} \\
        ${bam_output} \\
        > ${bam_format ? prefix + ".bam" : prefix + ".paf"}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def out    = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    """
    touch ${out}
    ${bam_format ? "touch ${prefix}.bam.bai" : ""}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.28
    END_VERSIONS
    """
}
