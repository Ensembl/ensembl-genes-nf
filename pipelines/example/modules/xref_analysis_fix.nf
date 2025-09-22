
/*
 * Problem: The object_xref table contains analysis_id foreign keys that reference
 * non-existent records in the analysis table, violating referential integrity.
 *
 * Solution: Update all orphaned cross-references to point to a default 'ensembl'
 * analysis entry, preserving the xref data while restoring database consistency.
 */

process FIX_XREFS {

    tag "${meta.id}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/xref_updates"

    input:
    tuple val(meta), val(core_db)

    output:
    tuple val(meta), path("${core_db}.metadata.log"), emit: logs

    script:
    """
    LOGFILE="${core_db}.metadata.log"
    SUMMARY="${core_db}.metadata.txt"

    echo "=== Fixing xrefs for ${meta.id} (${core_db}) ===" > "\$LOGFILE"
    echo "Server: ${meta.server}" >> "\$LOGFILE"
    echo "Core DB: ${core_db}" >> "\$LOGFILE"
    echo "Date: \$(date)" >> "\$LOGFILE"
    echo >> "\$LOGFILE"

    # Run queries individually, log result
    echo "Fixing xrefs in ${core_db} on ${meta.server}" >> "\$LOGFILE"
    ${meta.server} ${core_db} -e "SELECT COUNT(*) FROM object_xref WHERE analysis_id NOT IN (SELECT analysis_id FROM analysis);" >> "\$LOGFILE" 2>&1
    ${meta.server} ${core_db} -e "UPDATE object_xref SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = 'ensembl')WHERE analysis_id NOT IN (SELECT analysis_id FROM analysis);"  >> "\$LOGFILE" 2>&1
    echo "Xrefs fixed in ${core_db} on ${meta.server}" >> "\$LOGFILE"
    echo "=== Summary of changes ===" >> "\$LOGFILE"
    ${meta.server} ${core_db} -e "SELECT COUNT(*) FROM object_xref WHERE analysis_id NOT IN (SELECT analysis_id FROM analysis);" >> "\$LOGFILE" 2>&1
    echo "=== End of xref fix ===" >> "\$LOGFILE"
    """
}

