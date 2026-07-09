process MINIMAP2_WGA {
    tag "${meta.id}"
    label 'process_high'

    container "https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0"

    publishDir "${params.outdir}/${meta.id}/synteny", mode: 'copy', pattern: "*.paf"

    input:
    tuple val(meta), path(target_fasta)   // per-target assembly
    path reference                          // shared reference FASTA

    output:
    tuple val(meta), path("${prefix}.paf"), emit: paf
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Whole-genome alignment front-end (swappable). Must emit cs:Z long tags so the
    // hpp projection backend can lift coordinates exactly. Mirrors hpp run_minimap2:
    //   minimap2 -cx asm5 --cs=long -t <cpus> <ref> <target>
    def args   = task.ext.args   ?: '-cx asm5 --cs=long'
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        ${args} \\
        -t ${task.cpus} \\
        ${reference} \\
        ${target_fasta} \\
        > ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: 2.28-r1209
    END_VERSIONS
    """
}
