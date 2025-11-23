process BEDGRAPH_TO_BIGWIG {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::ucsc-bedgraphtobigwig=469"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:469--h9b8f530_0' :
        'biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0' }"

    publishDir "${params.outdir}/bigwigs", mode: 'copy'

    input:
    tuple val(meta), path(bedgraph)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Remove '.sorted.bedgraph' from the input filename to create the prefix
    def prefix = bedgraph.name.replaceFirst(/\.sorted\.bedgraph$/, '')
    """
    bedGraphToBigWig ${bedgraph} ${chrom_sizes} ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: \$(bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | sed 's/^.*v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    // Remove '.sorted.bedgraph' from the input filename to create the prefix
    def prefix = bedgraph instanceof List ? bedgraph[0].name.replaceFirst(/\.sorted\.bedgraph$/, '') : bedgraph.name.replaceFirst(/\.sorted\.bedgraph$/, '')
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: 469
    END_VERSIONS
    """
}
