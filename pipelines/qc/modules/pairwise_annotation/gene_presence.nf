process PAIRWISE_GENE_PRESENCE {
    tag { meta.id }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
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

        normalise_pairwise_tsv.py --input ${meta.id}_gene_presence.tsv --output ${meta.id}_gene_presence.normalised.tsv
        mv ${meta.id}_gene_presence.normalised.tsv ${meta.id}_gene_presence.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tgene_name\\tpresent_in_source_a\\tpresent_in_source_b\\tsource_a_gene_id\\tsource_b_gene_id\\n" > ${meta.id}_gene_presence.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
