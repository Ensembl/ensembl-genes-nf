process TRANSLON_DB {
    label 'process_ultra_high'

    container "community.wave.seqera.io/library/pip_duckdb_pandas_pyfaidx:2042eaa57c64430c"

    publishDir "${params.outdir}/translon_db", mode: 'copy'

    input:
    path input_root
    path genome_fasta
    path gencode_gtf
    path manifest

    output:
    path("translon_db/*"), emit: outputs

    script:
    def fasta_arg    = genome_fasta.name    != 'NO_FILE' ? "--fasta ${genome_fasta}"    : ""
    def gtf_arg      = gencode_gtf.name     != 'NO_FILE' ? "--gtf ${gencode_gtf}"       : ""
    def manifest_arg = manifest.name        != 'NO_FILE' ? "--manifest ${manifest}"     : ""
    """
    mkdir -p translon_db

    translon_db_standardise.py \\
        --input-root ${input_root} \\
        --out-dir translon_db \\
        ${fasta_arg} \\
        ${gtf_arg} \\
        ${manifest_arg}
    """

    stub:
    """
    mkdir -p translon_db
    echo '{"inputs":0,"translons":0}' > translon_db/translon_db_summary.json
    touch translon_db/translons.sqlite
    """
}
