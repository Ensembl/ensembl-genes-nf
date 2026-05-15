process PAIRWISE_GENE_PRESENCE {
    tag { meta.id }
    label 'process_low'

    container "python:3.11-slim"

    publishDir "${params.outdir}/qc/pairwise_annotation/${meta.id}/gene_presence",
        mode: 'copy',
        pattern: "*_gene_presence.tsv"

    input:
        tuple val(meta), path(source_a_gff), path(source_b_gff)
        val ensg_lookup

    output:
        tuple val(meta), path("${meta.id}_gene_presence.tsv"), emit: metrics
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def lookupArg = ensg_lookup ? "--ensg-lookup ${ensg_lookup}" : ''
        """
        set -euo pipefail

        compare_gene_presence.py \\
            --ensembl-gff ${source_a_gff} \\
            --cat-gff ${source_b_gff} \\
            --output ${meta.id}_gene_presence.tsv \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id} \\
            ${lookupArg}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tgene_name\\tpresent_in_ensembl\\tpresent_in_cat\\tensembl_gene_id\\tcat_gene_id\\n" > ${meta.id}_gene_presence.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
