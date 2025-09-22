
/*
 * Problem: Genes have canonical_transcript_id values pointing to deleted transcripts, 
 * and some genes exist without any transcripts (orphaned from projection cleanup).
 *
 * Solution: Delete orphaned genes with no transcripts, then update remaining broken 
 * canonical references to point to the longest coding transcript (or longest transcript).
 */

process ENSURE_CANONICAL_TX {
    tag "${meta.id}"

    // errorStrategy 'ignore'

    publishDir "${params.outdir}/canonical_fixes"

    input:
    tuple val(meta), val(core_db)
    tuple val(dummy), val(file) // Dummy input to ensure process runs after delete

    output:
    tuple val(meta), path("${core_db}.canonical_fix.txt"), emit: results
    tuple val(meta), path("${core_db}.canonical_fix.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.canonical_fix.log"
    SUMMARY="${core_db}.canonical_fix.txt"

    echo "=== Fixing canonical transcripts for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # Count broken canonicals before fix
    BROKEN_BEFORE=\$(${meta.server} ${core_db} -e "
    SELECT COUNT(*) FROM gene g
    WHERE g.canonical_transcript_id NOT IN (
        SELECT transcript_id FROM transcript WHERE gene_id = g.gene_id
    );" 2>> "\$LOGFILE" | tail -n 1)

    echo "Broken canonical transcripts before fix: \$BROKEN_BEFORE" >> "\$LOGFILE"

    if [ "\$BROKEN_BEFORE" -eq 0 ]; then
        echo "No broken canonical transcripts found" > "\$SUMMARY"
        exit 0
    fi

    # Identify and count orphaned genes (genes without transcripts)
    ORPHANED_GENES=\$(${meta.server} ${core_db} -e "
    SELECT COUNT(*) FROM gene g
    WHERE NOT EXISTS (SELECT 1 FROM transcript WHERE gene_id = g.gene_id);" 2>> "\$LOGFILE" | tail -n 1)

    echo "Orphaned genes (no transcripts): \$ORPHANED_GENES" >> "\$LOGFILE"

    # Log orphaned genes before deletion
    if [ "\$ORPHANED_GENES" -gt 0 ]; then
        echo "Logging orphaned genes..." >> "\$LOGFILE"
        ${meta.server} ${core_db} -e "
        SELECT CONCAT('ORPHANED: ', stable_id, ' (gene_id=', gene_id, ')') as deleted_gene
        FROM gene g
        WHERE NOT EXISTS (SELECT 1 FROM transcript WHERE gene_id = g.gene_id);" >> "\$LOGFILE" 2>&1
        
        # Delete orphaned genes
        echo "Deleting orphaned genes..." >> "\$LOGFILE"
        ${meta.server} ${core_db} -e "
        DELETE FROM gene 
        WHERE NOT EXISTS (SELECT 1 FROM transcript WHERE gene_id = gene.gene_id);" >> "\$LOGFILE" 2>&1
        
        echo "Deleted \$ORPHANED_GENES orphaned genes" >> "\$LOGFILE"
    fi

    # Now fix remaining broken canonical transcripts
    REMAINING_BROKEN=\$(${meta.server} ${core_db} -e "
    SELECT COUNT(*) FROM gene g
    WHERE g.canonical_transcript_id NOT IN (
        SELECT transcript_id FROM transcript WHERE gene_id = g.gene_id
    );" 2>> "\$LOGFILE" | tail -n 1)

    echo "Broken canonical transcripts after orphan cleanup: \$REMAINING_BROKEN" >> "\$LOGFILE"

    if [ "\$REMAINING_BROKEN" -gt 0 ]; then
        echo "Updating broken canonical transcripts..." >> "\$LOGFILE"
        ${meta.server} ${core_db} -e "
        UPDATE gene g
        SET g.canonical_transcript_id = (
            SELECT t.transcript_id
            FROM transcript t
            LEFT JOIN translation tl ON t.transcript_id = tl.transcript_id
            WHERE t.gene_id = g.gene_id
            ORDER BY 
                CASE WHEN tl.translation_id IS NOT NULL THEN 1 ELSE 2 END,
                COALESCE(tl.seq_end - tl.seq_start + 1, 0) DESC,
                (t.seq_region_end - t.seq_region_start + 1) DESC
            LIMIT 1
        )
        WHERE g.canonical_transcript_id NOT IN (
            SELECT transcript_id FROM transcript WHERE gene_id = g.gene_id
        );" >> "\$LOGFILE" 2>&1
    fi

    # Final verification
    BROKEN_AFTER=\$(${meta.server} ${core_db} -e "
    SELECT COUNT(*) FROM gene g
    WHERE g.canonical_transcript_id NOT IN (
        SELECT transcript_id FROM transcript WHERE gene_id = g.gene_id
    );" 2>> "\$LOGFILE" | tail -n 1)

    echo "Broken canonical transcripts after fix: \$BROKEN_AFTER" >> "\$LOGFILE"

    # Summary
    FIXED_COUNT=\$((\$REMAINING_BROKEN - \$BROKEN_AFTER))
    echo "Summary:" > "\$SUMMARY"
    echo "- Deleted \$ORPHANED_GENES orphaned genes" >> "\$SUMMARY"
    echo "- Fixed \$FIXED_COUNT canonical transcripts" >> "\$SUMMARY"
    
    if [ "\$BROKEN_AFTER" -eq 0 ]; then
        echo "- SUCCESS: All canonical transcripts now valid" >> "\$SUMMARY"
    else
        echo "- WARNING: \$BROKEN_AFTER canonical transcripts still broken" >> "\$SUMMARY"
    fi
    
    echo "Total genes processed: \$((\$BROKEN_BEFORE)) -> remaining issues: \$BROKEN_AFTER" >> "\$SUMMARY"
    """
}