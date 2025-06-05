
process GET_NCBI_TAXONOMY {
    errorStrategy 'ignore'

    publishDir "${params.outdir}/taxa", mode: 'copy'

    output:
    path("taxa.sqlite"), emit: taxa_db

    script:
    """
    get_ncbi_taxonomy.py -o taxa.sqlite
    """
}
