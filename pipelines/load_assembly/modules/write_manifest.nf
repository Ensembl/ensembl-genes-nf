process WRITE_MANIFEST {
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${outdir}", mode: 'copy', overwrite: true

    input:
    val  outdir
    path genome_fasta
    path synonyms_tsv
    path meta_json

    output:
    path 'output_manifest.json', emit: manifest

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    write_manifest.py \\
        --pipeline    load_assembly --version 1.0.0 \\
        --outdir      ${outdir} \\
        --genome_fasta ${genome_fasta} \\
        --synonyms_tsv ${synonyms_tsv} \\
        --meta_json    ${meta_json}
    """

    stub:
    """
    echo '{"pipeline":"load_assembly","version":"1.0.0","completed_at":"stub","outputs":[]}' > output_manifest.json
    """
}
