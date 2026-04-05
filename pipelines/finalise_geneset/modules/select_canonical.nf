// SELECT_CANONICAL
// For each gene, select the canonical transcript:
//   1. Longest CDS (sum of CDS feature lengths)
//   2. Tie-break: longest total transcript span
//   3. Tie-break: first in file order
// All transcripts are retained.  The canonical transcript gains the attribute
// canonical_transcript=1; all others receive canonical_transcript=0.

process SELECT_CANONICAL {
    label 'process_low'

    conda "conda-forge::python=3.11 bioconda::pybedtools=0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.10.0--py311h7b4e6c6_1' :
        'biocontainers/pybedtools:0.10.0--py311h7b4e6c6_1' }"

    publishDir path: "${params.outdir}/finalise_geneset", mode: 'copy', overwrite: true,
               pattern: '*.final.gff3'

    input:
    path  input_gff3

    output:
    path "*.final.gff3", emit: gff3
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.final.gff3'
    """
    select_canonical.py \\
        --gff3 ${input_gff3} \\
        --out ${out_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def out_name = input_gff3.name.replaceAll(/\.gff3$/, '') + '.final.gff3'
    """
    cp ${input_gff3} ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
