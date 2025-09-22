/*
 * Problem: Attribute and cross-reference tables contain records pointing to deleted genes, 
 * transcripts, and translations (likely poorly created or broken features from genome 
 * projection), creating orphaned metadata that wastes space and violates referential integrity.
 * 
 * Solution: Use LEFT JOIN queries to identify and DELETE all orphaned attributes, 
 * cross-references, synonyms, and unreferenced xrefs across the genomic feature hierarchy.
 */

process CLEAN_ORPHANED_METADATA {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/orphan_cleanup"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), path("${core_db}.orphan_cleanup.txt"), emit: results
    tuple val(meta), path("${core_db}.orphan_cleanup.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.orphan_cleanup.log"
    SUMMARY="${core_db}.orphan_cleanup.txt"

    echo "=== Orphaned Metadata Cleanup for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # Delete orphaned gene attributes
    ${meta.server} ${core_db} -e "DELETE ga FROM gene_attrib ga LEFT JOIN gene g ON ga.gene_id = g.gene_id WHERE g.gene_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Delete orphaned transcript attributes  
    ${meta.server} ${core_db} -e "DELETE ta FROM transcript_attrib ta LEFT JOIN transcript t ON ta.transcript_id = t.transcript_id WHERE t.transcript_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Delete orphaned translation attributes
    ${meta.server} ${core_db} -e "DELETE ta FROM translation_attrib ta LEFT JOIN translation tl ON ta.translation_id = tl.translation_id WHERE tl.translation_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Delete orphaned gene xrefs
    ${meta.server} ${core_db} -e "DELETE ox FROM object_xref ox LEFT JOIN gene g ON ox.ensembl_id = g.gene_id WHERE ox.ensembl_object_type = 'Gene' AND g.gene_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Delete orphaned transcript xrefs
    ${meta.server} ${core_db} -e "DELETE ox FROM object_xref ox LEFT JOIN transcript t ON ox.ensembl_id = t.transcript_id WHERE ox.ensembl_object_type = 'Transcript' AND t.transcript_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Delete orphaned translation xrefs
    ${meta.server} ${core_db} -e "DELETE ox FROM object_xref ox LEFT JOIN translation tl ON ox.ensembl_id = tl.translation_id WHERE ox.ensembl_object_type = 'Translation' AND tl.translation_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Clean up unreferenced external synonyms
    ${meta.server} ${core_db} -e "DELETE es FROM external_synonym es LEFT JOIN xref x ON es.xref_id = x.xref_id WHERE es.xref_id IS NOT NULL AND x.xref_id IS NULL;" >> "\$LOGFILE" 2>&1

    # Clean up unreferenced xrefs
    ${meta.server} ${core_db} -e "DELETE x FROM xref x WHERE NOT EXISTS (SELECT 1 FROM object_xref ox WHERE ox.xref_id = x.xref_id);" >> "\$LOGFILE" 2>&1

    echo "Cleanup completed: \$(date)" >> "\$LOGFILE"
    echo "Orphaned metadata cleanup completed for ${core_db}" > "\$SUMMARY"
    """
}