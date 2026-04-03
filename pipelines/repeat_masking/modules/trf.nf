// TRF — Tandem Repeat Finder
// Identifies tandem repeats in genomic sequences.
// Outputs a BED file of tandem repeat coordinates.

process TRF {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::trf=4.09.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trf:4.09.1--hec16e2b_3' :
        'biocontainers/trf:4.09.1--hec16e2b_3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.trf.bed"), emit: bed
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    // TRF params: match, mismatch, delta, PM, PI, minscore, maxperiod
    def args   = task.ext.args ?: '2 7 7 80 10 50 500'
    """
    trf ${fasta} ${args} -f -d -m -ngs > ${prefix}.trf.dat

    trf_to_bed.py --dat ${prefix}.trf.dat --out ${prefix}.trf.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trf: \$(trf 2>&1 | head -1 | awk '{print \$2}' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.trf.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trf: 4.09.1
    END_VERSIONS
    """
}
