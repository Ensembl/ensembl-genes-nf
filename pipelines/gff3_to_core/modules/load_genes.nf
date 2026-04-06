process LOAD_GFF3_TO_CORE {
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pymysql=1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11--h2ad013b_0_cp311' :
        'biocontainers/python:3.11--h2ad013b_0_cp311' }"

    publishDir path: "${params.outdir}/gff3_to_core", mode: 'copy', overwrite: true

    input:
    path gff3
    path genome_fai
    path synonyms_tsv

    output:
    path 'load_stats.json', emit: stats
    path 'versions.yml',    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fai_arg     = (genome_fai.name    != 'NO_FILE_FAI') ? "--genome-fai ${genome_fai}"     : ''
    def syn_arg     = (synonyms_tsv.name  != 'NO_FILE_SYN') ? "--synonyms-tsv ${synonyms_tsv}" : ''
    """
    load_gff3_to_core.py \\
        --gff3          ${gff3} \\
        --host          ${params.db_host} \\
        --port          ${params.db_port} \\
        --user          ${params.db_user} \\
        --password      ${params.db_password} \\
        --dbname        ${params.db_name} \\
        --analysis      ${params.analysis_logic_name} \\
        --coord-system  ${params.coord_system} \\
        --assembly      ${params.assembly} \\
        --species-id    ${params.species_id} \\
        --species-name  ${params.species_name} \\
        ${fai_arg} \\
        ${syn_arg} \\
        --stats-json    load_stats.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pymysql: \$(python3 -c "import pymysql; print(pymysql.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{"genes":10,"transcripts":25,"coding_transcripts":20,"exons":80,"translations":20,"skipped_genes":0}' \\
        > load_stats.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        pymysql: 1.1.0
    END_VERSIONS
    """
}
