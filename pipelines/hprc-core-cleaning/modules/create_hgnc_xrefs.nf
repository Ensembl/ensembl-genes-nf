
/*
 * Problem: Genes with HGNC identifiers in their descriptions lack corresponding 
 * xref records, preventing proper display_xref_id assignment.
 *
 * Solution: Parse HGNC accessions and display labels from gene descriptions to 
 * create missing xref entries in the external database (external_db_id 1100).
 */

process CREATE_HGNC_XREFS {
    
    tag "${meta.id}"
    
    publishDir "${params.outdir}/gene_display_xref_updates"
    
    input:
    tuple val(meta), val(core_db)
    
    output:
    tuple val(meta), val(core_db), path("${core_db}.xref_creation.log"), emit: created_xrefs
    
    script:
    """
    LOGFILE="${core_db}.xref_creation.log"
    
    echo "=== Creating HGNC xrefs for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"
    
    ${meta.server} ${core_db} -e "
    INSERT IGNORE INTO xref (external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text)
    SELECT 
        1100 as external_db_id,
        CONCAT('HGNC:', SUBSTRING(description, 
            LOCATE('HGNC:', description) + 5, 
            LOCATE(']', description, LOCATE('HGNC:', description)) - LOCATE('HGNC:', description) - 5
        )) as dbprimary_acc,
        SUBSTRING(description,
            LOCATE('parent_gene_display_xref=', description) + 25,
            LOCATE(';', description, LOCATE('parent_gene_display_xref=', description)) - LOCATE('parent_gene_display_xref=', description) - 25
        ) as display_label,
        0 as version,
        SUBSTRING(description, 1, LOCATE('[', description) - 1) as description,
        'DIRECT' as info_type,
        'Generated via ensembl_manual' as info_text
    FROM gene 
    WHERE display_xref_id IS NULL 
      AND description LIKE '%Source:HGNC Symbol;Acc:HGNC:%'
      AND description LIKE '%parent_gene_display_xref=%';
    " >> "\$LOGFILE" 2>&1
    
    echo "Xref creation exit code: \$?" >> "\$LOGFILE"
    """
}
