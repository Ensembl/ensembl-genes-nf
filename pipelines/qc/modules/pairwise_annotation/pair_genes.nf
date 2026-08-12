process PAIR_GENES {
    tag { meta.id }
    label 'process_medium'

    container "https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1"

    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*.tsv"
    publishDir "${params.outdir}/qc/pairwise_annotation",
        mode: 'copy',
        pattern: "*.log"

    input:
        tuple val(meta), path(source_a_gff, stageAs: 'source_a.gff3'), path(source_b_gff, stageAs: 'source_b.gff3'), path(assembly_report)

    output:
        tuple val(meta), path("${meta.id}.gene_pairs_rbh.tsv"), emit: rbh
        tuple val(meta), path("${meta.id}.gene_pairs_all.tsv"), emit: all_pairs
        tuple val(meta), path("${meta.id}.assembly_summary.tsv"), emit: summary
        tuple val(meta), path("${meta.id}.unmatched_genes.tsv"), optional: true, emit: unmatched
        tuple val(meta), path("${meta.id}.biotype_stats.tsv"), optional: true, emit: biotype_stats
        tuple val(meta), path("${meta.id}.log"), emit: log
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def normalization = params.pairwise_contig_normalization ?: 'none'
        """
        set -euo pipefail

        mkdir -p dummy_dirs/source_a dummy_dirs/source_b
        ln -s "\$(readlink -f source_a.gff3)" dummy_dirs/source_a/
        ln -s "\$(readlink -f source_b.gff3)" dummy_dirs/source_b/

        printf "Assembly Accession,Sample Name\\n%s,%s\\n" "${meta.id}" "${meta.id}" > assemblies.csv
        printf "Sample Name,GFF File\\n%s,%s\\n" "${meta.id}" "source_b.gff3" > source_b_index.csv

        hprc_ensembl_cat_overlap.py \\
            --ensembl-gff source_a.gff3 \\
            --cat-gff source_b.gff3 \\
            --ensembl-gff-dir dummy_dirs/source_a \\
            --cat-anno-dir dummy_dirs/source_b \\
            --cat-index source_b_index.csv \\
            --assemblies-index assemblies.csv \\
            --output-prefix ${meta.id} \\
            --contig-normalization ${normalization} \\
            ${args} \\
            2>&1 | tee ${meta.id}.log

        for table in ${meta.id}.gene_pairs_all.tsv ${meta.id}.gene_pairs_rbh.tsv ${meta.id}.assembly_summary.tsv; do
            normalise_pairwise_tsv.py --input "\$table" --output "\$table.normalised"
            mv "\$table.normalised" "\$table"
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """

    stub:
        """
        printf "assembly_accession\\tsample_name\\tsource_a_gene_id\\tsource_a_biotype\\tsource_b_gene_id\\tsource_b_biotype\\toverlap_bp\\tfrac_source_a_covered\\tfrac_source_b_covered\\tclassification\\tclassification_detailed\\tis_rbh\\n" > ${meta.id}.gene_pairs_all.tsv
        printf "assembly_accession\\tsample_name\\tsource_a_gene_id\\tsource_a_biotype\\tsource_b_gene_id\\tsource_b_biotype\\toverlap_bp\\tfrac_source_a_covered\\tfrac_source_b_covered\\tclassification\\tclassification_detailed\\tis_rbh\\n" > ${meta.id}.gene_pairs_rbh.tsv
        printf "assembly_accession\\tsample_name\\tn_source_a_genes\\tn_source_b_genes\\tn_rbh_pairs\\n" > ${meta.id}.assembly_summary.tsv
        touch ${meta.id}.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.11.0
        END_VERSIONS
        """
}
