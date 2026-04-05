// MERGE_REPEATS
// Merge BED files from RepeatMasker, RED, TRF, and DUST into a single
// non-redundant sorted GFF3 + BED file ready for downstream use.

process MERGE_REPEATS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/repeats", mode: 'copy', pattern: "*.repeats.gff3"

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h18e979d_0' :
        'biocontainers/pybedtools:0.10.0--py311h18e979d_0' }"

    input:
    tuple val(meta), path(beds)   // all BED files collected into one list

    output:
    tuple val(meta), path("*.repeats.gff3"), emit: gff3
    tuple val(meta), path("*.repeats.bed"),  emit: bed
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args   = task.ext.args   ?: ''
    """
    merge_repeats.py \\
        --beds     ${beds} \\
        --out-gff3 ${prefix}.repeats.gff3 \\
        --out-bed  ${prefix}.repeats.bed \\
        --sample   ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)" 2>/dev/null || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.repeats.gff3 ${prefix}.repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
        pybedtools: 0.10.0
    END_VERSIONS
    """
}
