// modules/local/audit_run.nf
process AUDIT_RUN {
    tag "$db_name"

    publishDir {
        "${params.output_dir}/${db_name}/audit"
    }, mode: 'copy',
    saveAs: { name -> name.startsWith("${db_name}.") ? name : null }

    input:
        tuple val(db_name),
              path(ref_gff),
              path(target_gff),
              val(mapping_session_id),
              path(locus_comparison),
              path(decisions_tsv),
              path(score_evidence_tsv),
              path(executable_sql),
              path(dry_run_sql)

    output:
    	tuple path("${db_name}.stable_id_audit.txt"),
    	      path("${db_name}.missing_genes.tsv"),
    	      path("${db_name}.coordinate_mapped_genes.tsv"),
    	      path("${db_name}.new_genes.tsv"),
    	      emit: audit

    script:
    """
    python3 ${projectDir}/bin/audit_stable_id_run.py \
        --decisions-tsv ${decisions_tsv} \
        --ref-gff ${ref_gff} \
        --target-gff ${target_gff} \
        --gene-locus-comparison ${locus_comparison} \
        --output-dir audit_tables \
        --limit ${params.audit_limit} \
        > ${db_name}.stable_id_audit.txt

    cp audit_tables/missing_genes.tsv ${db_name}.missing_genes.tsv
    cp audit_tables/coordinate_mapped_genes.tsv ${db_name}.coordinate_mapped_genes.tsv
    cp audit_tables/new_genes.tsv ${db_name}.new_genes.tsv
    """
}
