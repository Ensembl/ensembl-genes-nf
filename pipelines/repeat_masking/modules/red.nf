// RED — Rapid Repeat Element Detector
// De novo repeat masking without a library. Fast but less sensitive than RepeatMasker.
// Outputs a BED file of repeat coordinates.

process RED {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "bioconda::red=05.2022.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/red:05.2022.02--h9ee0642_1' :
        'biocontainers/red:05.2022.02--h9ee0642_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.rpt.bed"), emit: bed
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    mkdir -p red_out

    Red \\
        -gnm ${fasta} \\
        -dir red_out \\
        -msk red_out \\
        ${args}

    red_to_bed.py --rpt red_out/*.rpt --out ${prefix}.rpt.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        red: \$(Red 2>&1 | head -1 | sed 's/Red v//' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.rpt.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        red: 05.2022.02
    END_VERSIONS
    """
}
