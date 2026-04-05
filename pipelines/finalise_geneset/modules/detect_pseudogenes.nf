// DETECT_PSEUDOGENES
// Flag pseudogenes by two criteria:
//   1. Single-exon coding genes where >max_repeat_cds_coverage of the CDS is
//      covered by repeat features → biotype 'processed_pseudogene'
//   2. Multi-exon coding genes where ALL introns are frameshifted (very short)
//      → biotype 'pseudogene'

process DETECT_PSEUDOGENES {
    label 'process_medium'

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h7b4e6c6_1' :
        'biocontainers/pybedtools:0.10.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/finalise_geneset", mode: 'copy', overwrite: true,
               pattern: '*.pseudogenes.gff3'

    input:
    path  filtered_gff3
    path  repeat_gff3
    val   max_repeat_cds_coverage

    output:
    path "*.pseudogenes.gff3", emit: gff3
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = filtered_gff3.name.replaceAll(/\.gff3$/, '') + '.pseudogenes.gff3'
    """
    detect_pseudogenes.py \\
        --gff3 ${filtered_gff3} \\
        --repeats ${repeat_gff3} \\
        --out ${out_name} \\
        --max-repeat-coverage ${max_repeat_cds_coverage} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = filtered_gff3.name.replaceAll(/\.gff3$/, '') + '.pseudogenes.gff3'
    """
    cp ${filtered_gff3} ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
