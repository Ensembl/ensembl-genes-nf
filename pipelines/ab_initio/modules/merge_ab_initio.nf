// MERGE_AB_INITIO
// Concatenate per-chunk ab initio GFF3 files and renumber IDs sequentially.

process MERGE_AB_INITIO {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${params.outdir}/ab_initio", mode: 'copy', overwrite: true,
               pattern: '*.merged.gff3'

    input:
    path gff3_files

    output:
    path "ab_initio.merged.gff3", emit: gff3
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    merge_ab_initio_gff3.py \\
        --inputs ${gff3_files} \\
        --out    ab_initio.merged.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '##gff-version 3\\n' > ab_initio.merged.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
    END_VERSIONS
    """
}
