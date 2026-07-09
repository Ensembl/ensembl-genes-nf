process BUILD_SYNTENY {
    tag "${meta.id}"
    label 'process_medium'

    container "${params.hpp_container}"

    input:
    tuple val(meta), path(paf), path(target_fasta)   // per-target: WGA PAF + assembly
    path reference                                     // shared reference FASTA
    path annotation                                    // shared reference GFF3

    output:
    tuple val(meta), path("${prefix}_state"), emit: state
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Parse the external PAF into syntenic blocks, fill synteny gaps, detect sex
    // chromosomes. The main whole-genome alignment is NOT run here (--paf supplies it).
    def args   = task.ext.args   ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    hpp build-synteny \\
        --ref-fasta ${reference} \\
        --ref-gff ${annotation} \\
        --target-fasta ${target_fasta} \\
        --paf ${paf} \\
        --output-gff ${prefix}.placeholder.gff3 \\
        --output-stats ${prefix}.placeholder.stats.json \\
        --threads ${task.cpus} \\
        ${args} \\
        --out-state ${prefix}_state

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: \$(hpp --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_state
    touch ${prefix}_state/syntenic_blocks.jsonl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpp: 0.1.0
    END_VERSIONS
    """
}
