// DETECT_READTHROUGH
// Flag transcripts that span two independent gene loci.
// A readthrough transcript's CDS overlaps with coding exons of two separate
// genes on the same strand.  Sets biotype to 'readthrough_transcript' on the
// transcript, and 'readthrough' on the gene if ALL its transcripts are flagged.

process DETECT_READTHROUGH {
    label 'process_low'

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h7b4e6c6_1' :
        'biocontainers/pybedtools:0.10.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/finalise_geneset", mode: 'copy', overwrite: true,
               pattern: '*.readthrough.gff3'

    input:
    path  input_gff3
    val   max_gap

    output:
    path "*.readthrough.gff3", emit: gff3
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.readthrough.gff3'
    """
    detect_readthrough.py \\
        --gff3 ${input_gff3} \\
        --out ${out_name} \\
        --max-gap ${max_gap} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.readthrough.gff3'
    """
    cp ${input_gff3} ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
