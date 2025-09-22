
/*
 * Problem: Genes have NULL display_xref_id values despite having HGNC identifiers 
 * and display labels embedded in their description fields from genome projection.
 *
 * Solution: Parse HGNC accessions and display labels from gene descriptions to 
 * match existing xref records and update the display_xref_id foreign key.
 */

process UPDATE_GENE_DISPLAY_XREF_IDS {
    
    tag "${meta.id}"
    
    publishDir "${params.outdir}/gene_display_xref_updates"
    
    input:
    tuple val(meta), val(core_db), path(xref_log)
    
    output:
    tuple val(meta), val(core_db), path("${core_db}.gene_update.log"), emit: updated_genes
    
    script:
    """
    LOGFILE="${core_db}.gene_update.log"
    
    echo "=== Updating gene display_xref_id for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"
    
    ${meta.server} ${core_db} -e "
    UPDATE gene g
    JOIN xref x ON (
        x.dbprimary_acc = CONCAT('HGNC:', SUBSTRING(g.description, 
            LOCATE('HGNC:', g.description) + 5, 
            LOCATE(']', g.description, LOCATE('HGNC:', g.description)) - LOCATE('HGNC:', g.description) - 5
        ))
        AND x.display_label = SUBSTRING(g.description,
            LOCATE('parent_gene_display_xref=', g.description) + 25
        )
        AND x.external_db_id = 1100
    )
    SET g.display_xref_id = x.xref_id
    WHERE g.display_xref_id IS NULL 
      AND g.description LIKE '%Source:HGNC Symbol;Acc:HGNC:%'
      AND g.description LIKE '%parent_gene_display_xref=%';
    " >> "\$LOGFILE" 2>&1
    
    echo "Gene update exit code: \$?" >> "\$LOGFILE"
    """
}