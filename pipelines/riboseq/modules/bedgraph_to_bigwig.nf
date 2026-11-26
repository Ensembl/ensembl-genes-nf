process BEDGRAPH_TO_BIGWIG {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::ucsc-bedgraphtobigwig=469"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:469--h9b8f530_0' :
        'biocontainers/ucsc-bedgraphtobigwig:469--h9b8f530_0' }"

    input:
    tuple val(meta), path(bedgraph)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bw"), optional: true, emit: bigwig
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Handle both single and multiple bedgraph files (stranded output)
    if (bedgraph instanceof List) {
        // Stranded output: process each bedgraph file
        def commands = bedgraph.collect { bg ->
            def prefix = bg.name.replaceFirst(/\.sorted\.bedgraph$/, '')
            "if [ -s '${bg}' ]; then bedGraphToBigWig '${bg}' '${chrom_sizes}' '${prefix}.bw'; else echo 'Skipping empty bedGraph: ${bg}' >&2; fi"
        }.join('\n    ')
        """
        ${commands}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedGraphToBigWig: \$(bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | sed 's/^.*v//; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        // Single bedgraph file
        def prefix = bedgraph.name.replaceFirst(/\.sorted\.bedgraph$/, '')
        """
        if [ -s '${bedgraph}' ]; then
            bedGraphToBigWig '${bedgraph}' '${chrom_sizes}' '${prefix}.bw'
        else
            echo 'Skipping empty bedGraph: ${bedgraph}' >&2
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedGraphToBigWig: \$(bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | sed 's/^.*v//; s/ .*\$//')
        END_VERSIONS
        """
    }

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
