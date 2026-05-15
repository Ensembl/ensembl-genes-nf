process CANONICAL_ORF_DB {
    label 'process_high'

    container "oras://community.wave.seqera.io/library/pip_duckdb_pandas:57c3741f55d53490"

    publishDir "${params.outdir}/canonical_orf_db", mode: 'copy'

    input:
    path input_root
    val genome_fasta
    val gencode_gtf

    output:
    path("canonical_orf_db/*"), emit: outputs

    script:
    def fasta_arg = genome_fasta ? "--fasta ${genome_fasta}" : ""
    def gtf_arg = gencode_gtf ? "--gtf ${gencode_gtf}" : ""
    """
    mkdir -p canonical_orf_db

    canonical_orf_standardise.py \\
        --input-root ${input_root} \\
        --out-dir canonical_orf_db \\
        ${fasta_arg} \\
        ${gtf_arg}
    """

    stub:
    """
    mkdir -p canonical_orf_db
    echo '{"inputs":0,"canonical_orfs":0}' > canonical_orf_db/canonical_orf_db_summary.json
    touch canonical_orf_db/canonical_orfs.sqlite
    """
}
