process FETCH_REFSEQ {

    tag "${meta.id}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta)

    output:
    tuple val(meta),
          path("refseq_data/**/*_ensembl.gff3"),
          path("refseq_data/**/*_genomic_ensembl.fna"),
          path("refseq_data/**/*_assembly_report.txt"),
          emit: refseq

    script:
    def script = "${ensembl_genes_repo}/src/python/ensembl/genes/ensembl_loading/gff_cli.py"

    """
    python ${script} \
        gff-loader \
        --log-file gff-loader.log \
        refseq \
        run \
        --base-dir refseq_data \
        --assembly-acc ${meta.id}
    """
}
