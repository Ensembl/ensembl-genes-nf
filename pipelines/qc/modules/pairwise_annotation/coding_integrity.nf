process PAIRWISE_CODING_INTEGRITY {
    tag { meta.id }
    label 'process_low'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*_coding_integrity.tsv"

    input:
        tuple val(meta), path(source_a_gff), path(source_b_gff), path(rbh_pairs)

    output:
        tuple val(meta), path("${meta.id}_coding_integrity.tsv"), emit: metrics
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail

        assess_coding_integrity.py \\
            --ensembl-gff ${source_a_gff} \\
            --cat-gff ${source_b_gff} \\
            --rbh-pairs ${rbh_pairs} \\
            --output ${meta.id}_coding_integrity.tsv \\
            --assembly-accession ${meta.id} \\
            --sample-name ${meta.id}

        normalise_pairwise_tsv.py --input ${meta.id}_coding_integrity.tsv --output ${meta.id}_coding_integrity.normalised.tsv
        mv ${meta.id}_coding_integrity.normalised.tsv ${meta.id}_coding_integrity.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tsource_a_gene_id\\tsource_b_gene_id\\tclassification\\n" > ${meta.id}_coding_integrity.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
