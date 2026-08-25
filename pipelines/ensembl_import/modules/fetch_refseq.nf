process FETCH_REFSEQ {
    label 'process_medium'

    tag "${meta.id}"

    publishDir { "${params.outdir}/refseq/${meta.id}" }, mode: 'copy', overwrite: true

    input:
    val meta

    output:
    tuple val(meta),
          path("refseq_data/**/*_ensembl.gff3"),
          path("refseq_data/**/*_genomic_ensembl.fna"),
          path("refseq_data/**/*_assembly_report.txt"),
          emit: refseq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gff-loader \
        --log-file gff-loader.log \
        refseq \
        run \
        ${args} \
        --base-dir refseq_data \
        --assembly-acc ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | awk '{print \$2}')
        ensembl-genes: \$(python -c 'from importlib.metadata import version; print(version("ensembl-genes"))')

    END_VERSIONS
    """

    stub:
    """
    mkdir -p refseq_data/${meta.id}
    touch refseq_data/${meta.id}/${meta.id}_ensembl.gff3
    touch refseq_data/${meta.id}/${meta.id}_genomic_ensembl.fna
    touch refseq_data/${meta.id}/${meta.id}_assembly_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
        ensembl-genes: stub
    END_VERSIONS
    """
}
