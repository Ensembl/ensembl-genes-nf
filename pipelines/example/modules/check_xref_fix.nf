
/*
 * Problem: After creating xrefs and updating display_xref_id values, need to verify 
 * that the fix worked and genes now have proper display xref assignments.
 *
 * Solution: Count genes with/without display_xref_id and show sample assignments 
 * to confirm the HGNC xref creation and linking process succeeded.
 */

process VERIFY_XREF_FIX {
    
    tag "${meta.id}"
    
    publishDir "${params.outdir}/gene_display_xref_verification"
    
    input:
    tuple val(meta), val(core_db), path(update_log)
    
    output:
    tuple val(meta), path("${core_db}.verification.log"), emit: verification
    
    script:
    """
    LOGFILE="${core_db}.verification.log"
    
    echo "=== Verification of xref fix for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"
    
    ${meta.server} ${core_db} -e "SELECT COUNT(*) as genes_still_without_display_xref FROM gene WHERE display_xref_id IS NULL;" >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "SELECT COUNT(*) as genes_now_with_display_xref FROM gene WHERE display_xref_id IS NOT NULL;" >> "\$LOGFILE" 2>&1
    
    ${meta.server} ${core_db} -e "
    SELECT g.stable_id, x.dbprimary_acc, x.display_label
    FROM gene g 
    JOIN xref x ON g.display_xref_id = x.xref_id 
    WHERE x.info_text = 'Generated via ensembl_manual'
    LIMIT 5;
    " >> "\$LOGFILE" 2>&1
    """
}