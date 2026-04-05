// FLAG_SELENOPROTEINS
// Match genes against a selenoprotein database.  When the selenoprotein FASTA
// is provided (i.e. the input is not the NO_FILE sentinel), genes whose Name
// attribute matches a selenoprotein sequence ID have their biotype set to
// 'selenoprotein'.  When NO_FILE is passed the input GFF3 is copied unchanged.

process FLAG_SELENOPROTEINS {
    label 'process_medium'

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h7b4e6c6_1' :
        'biocontainers/pybedtools:0.10.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/finalise_geneset", mode: 'copy', overwrite: true,
               pattern: '*.selenoproteins.gff3'

    input:
    path  input_gff3
    path  selenoprotein_fasta

    output:
    path "*.selenoproteins.gff3", emit: gff3
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.selenoproteins.gff3'
    """
    flag_selenoproteins.py \\
        --gff3 ${input_gff3} \\
        --proteins ${selenoprotein_fasta} \\
        --out ${out_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.selenoproteins.gff3'
    """
    cp ${input_gff3} ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
