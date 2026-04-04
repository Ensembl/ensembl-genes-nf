process WRITE_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${outdir}", mode: 'copy', overwrite: true

    input:
    val  outdir
    path gff3

    output:
    path 'output_manifest.json', emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    write_manifest.py \\
        --pipeline ab_initio --version 1.0.0 \\
        --outdir ${outdir} --gff3 ${gff3}
    """

    stub:
    """
    echo '{"pipeline":"ab_initio","version":"1.0.0","completed_at":"stub","outputs":[]}' > output_manifest.json
    """
}
