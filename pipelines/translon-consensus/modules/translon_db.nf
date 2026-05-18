process TRANSLON_DB {
    label 'process_ultra_high'

    container "community.wave.seqera.io/library/pip_duckdb_pandas_pyfaidx:2042eaa57c64430c"

    publishDir "${params.outdir}/translon_db", mode: 'copy'

    input:
    path input_root
    path genome_fasta
    path gencode_gtf

    output:
    path("translon_db/*"), emit: outputs

    script:
    def fasta_arg = genome_fasta ? "--fasta ${genome_fasta}" : ""
    def gtf_arg = gencode_gtf ? "--gtf ${gencode_gtf}" : ""
    """
    mkdir -p translon_db

    translon_db_standardise.py \\
        --input-root ${input_root} \\
        --out-dir translon_db \\
        ${fasta_arg} \\
        ${gtf_arg}
    """

    stub:
    """
    mkdir -p translon_db
    echo '{"inputs":0,"translons":0}' > translon_db/translon_db_summary.json
    touch translon_db/translons.sqlite
    """
}
